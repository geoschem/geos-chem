! $Id: tpcore_mod.f,v 1.1 2009/09/16 14:06:00 bmy Exp $
      MODULE TPCORE_MOD
!
!******************************************************************************
!  Module TPCORE_MOD contains the TPCORE transport subroutine package by
!  S-J Lin, version 7.1. (bmy, 7/16/01, 9/18/07)
!
!  Module Routines:
!  ============================================================================
!  (1 ) TPCORE    : TPCORE driver program
!  (2 ) COSA      : TPCORE intermediate subroutine
!  (3 ) COSC      : TPCORE intermediate subroutine
!  (4 ) FCT3D     : TPCORE intermediate subroutine
!  (5 ) FILEW     : TPCORE intermediate subroutine
!  (6 ) FILNS     : TPCORE intermediate subroutine
!  (7 ) FXPPM     : TPCORE intermediate subroutine
!  (8 ) FYPPM     : TPCORE intermediate subroutine
!  (9 ) FZPPM     : TPCORE intermediate subroutine
!  (10) HILO      : TPCORE intermediate subroutine
!  (11) HILO3D    : TPCORE intermediate subroutine
!  (12) LMTPPM    : TPCORE intermediate subroutine
!  (13) QCKXYZ    : TPCORE intermediate subroutine
!  (14) XADV      : TPCORE intermediate subroutine
!  (15) XMIST     : TPCORE intermediate subroutine
!  (16) XTP       : TPCORE intermediate subroutine
!  (17) YMIST     : TPCORE intermediate subroutine
!  (18) YTP       : TPCORE intermediate subroutine
!  (19) PRESS_FIX : Wrapper for pressure-fixer subroutine DYN0
!  (20) DYN0      : Implements pressure fix for mass fluxes in TPCORE
!  (21) PFILTR    : Applies pressure filter to ALFA and BETA mass fluxes 
!  (22) LOCFLT    : Local pressure filter -- called from PFILTR
!  (23) POLFLT    : Polar pressure filter -- called from PFILTR
!  (24) DIAG_FLUX : Computes TPCORE mass fluxes for ND24, ND25, ND26 diags
!  (25) MASSFLUX  : Compute size-resolved aerosol mass flux based on
!                   fluxes of aerosol number
!
!  Functions:
!  ============================================================================
!  (1 ) IDFLUX    : Take tracer ID as input and return a switch determining
!                   if the flux shall be calculated for aerosol mass
!  (2 ) NBIN      : For size-resolved aerosol, take tracer number and return
!                   a size-bin number
!
!  GEOS-CHEM modules referenced by tagged_co_mod.f
!  ============================================================================
!  (1 ) diag_mod.f       : Module containing GEOS-CHEM diagnostic arrays
!  (2 ) dao_mod.f        : Module containing DAO met field arrays
!  (3 ) global_ch4_mod.f : Module containing routines to read 3-D CH4 field
!  (4 ) grid_mod.f       : Module containing horizontal grid information
!  (5 ) pressure_mod.f   : Module containing routines to compute P(I,J,L)
!  (6 ) time_mod.f       : Module containing routines to compute date & time
!
!  NOTES:
!  (1 ) The TPCORE subroutines have not been modified, except to replace
!        obsolete parallel loop directives.  It is more convenient to place
!        all of the TPCORE subroutines into a single module, this reduces
!        clutter. (bmy, 7/16/01)
!  (2 ) All parallel loops are now specified with OpenMP directives, 
!        for cross-platform compatibility. (bmy, 7/16/01)
!  (3 ) The routines in TPCORE_MOD have been validated against the previous
!        version (Code_4.16). (bmy, 7/16/01)
!  (4 ) Updated comments (bmy, 9/4/01)
!  (5 ) Removed obsolete code from 7/12/01.  Also implemented pressure-fix
!        subroutines PRESS_FIX, DYN0, PFLITR, LOCFLT, POLFLT.  (bmy, 10/9/01)
!  (6 ) Now use PSC2 instead of PS in subroutine DYN0.  Also delineate the
!        first-time header text with horizontal lines. (bdf, bmy, 4/15/02)
!  (7 ) Now zero XMASS_PF and YMASS_PF arrays on every call to TPCORE.  
!        This will avoid floating-point exceptions on the Alpha platform.
!        (bmy, 4/18/02)
!  (8 ) Now divide module header into MODULE PRIVATE, MODULE VARIABLES, and
!        MODULE ROUTINES sections.  Updated comments (bmy, 5/28/02)
!  (9 ) Deleted obsolete code from 4/02.  (bdf, bmy, 8/22/02)
!  (10) Minor bug fix for ALPHA platform: delete extra comma in format 
!        statement 2 in routine TPCORE.  Bug fix: now stop the run if NDT is
!        too large.  This makes sure we don't violate the Courant limit.
!        (bmy, 11/22/02)
!  (11) Also add output for the SUN/Sparc platform.  Rename DEC_COMPAQ to
!        COMPAQ.  Also assume that all platforms other than CRAY use OPENMP 
!        parallelization commands (bmy, 3/23/03)
!  (12) Now references "grid_mod.f" and "time_mod.f" (bmy, 3/24/03)
!  (13) Now print output for IBM/AIX platform in "tpcore" (gcc, bmy, 6/27/03)
!  (14) Remove obsolete code for CO-OH parameterization (bmy, 6/24/05)
!  (15) Bug fix in DIAG_FLUX: now dimension FX, FX properly (bmy, 7/21/05)
!  (16) Now print output for IFORT compiler in "tpcore" (bmy, 10/18/05)
!  (17) Remove support for LINUX_IFC & LINUX_EFC compilers (bmy, 8/4/06)
!  (18) Corrected mass flux diagnostics (phs, 9/18/07)
!  (19) Add subroutine MASSFLUX, function IDFLUX and NBIN and modify TPCORE
!       for aerosol number and mass advection by adding arrays for saving
!       fluxes and tracer mixing ratio fields arrays. (win, 7/17/09)
!  (19) Modify subroutine YTP for a complete saving of diagnostic array fy1_tp 
!      (win,7/17/09)
!  (20) Modify subroutine FZPPM extending the loop for fz1_tp (win, 7/17/09)
!  01 Mar 2012 - R. Yantosca - Now reference new grid_mod.F90
!******************************************************************************
!
      !=================================================================
      ! MODULE PRIVATE DECLARATIONS -- keep certain internal variables 
      ! and routines from being seen outside "tpcore_mod.f"
      !=================================================================

      USE ERROR_MOD, ONLY : ALLOC_ERR

      ! Make everything PRIVATE ...
      PRIVATE

      ! ... except this routine
      PUBLIC :: TPCORE

      !=================================================================
      ! MODULE ROUTINES -- follow below the "CONTAINS" statement 
      !=================================================================
      CONTAINS

!------------------------------------------------------------------------------

! ****6***0*********0*********0*********0*********0*********0**********72
      subroutine tpcore( IGD,  Q,    PS1,   PS2,   U,     V,
     &                   W,    NDT,  IORD,  JORD,  KORD,  NC,
     &                   IM,   JM,   j1,    NL,    AP,    BP,   
     &                   PT,   AE,   FILL,  MFCT,  Umax,  State_Met )
!****6***0*********0*********0*********0*********0*********0**********72
 
! TransPort module for Goddard Chemistry Transport Model (G-CTM), Goddard
! Earth Observing System General Circulation Model (GEOS-GCM), and Data
! Assimilation System (GEOS-DAS).
 
! Purpose: perform the transport of  3-D mixing ratio fields using 
!          externally specified winds on the hybrid Eta-coordinate.
!          One call to tpcore updates the 3-D mixing ratio
!          fields for one time step (NDT). [vertical mass flux is computed
!          internally using a center differenced hydrostatic mass
!          continuity equation].
 
! Schemes: Multi-dimensional Flux Form Semi-Lagrangian (FFSL) schemes
!          (Lin and Rood 1996, MWR) with a modified MFCT option (Zalesak 1979).
 
! Multitasking version: 7.1
! Last modified: Sept 2, 1999
! Changes from version 7.m: large-time-step bug in xtp fixed.
! Suggested compiler options:
! CRAY f77 compiler:  cf77 -Zp -c -Wd'-dec' -Wf' -a stack -exm'
! CRAY f90 compiler:  f90 -c -eZ -DCRAY -Dmultitask
! SGI Origin: f77 -c -DSGI -Dmultitask -r8 -64 -O3 -mips4 -mp
!             loader: f77 -64 -mp
!
! Send comments/suggestions to
!
!                 S.-J. Lin
! Address:
!                 Code 910.3, NASA/GSFC, Greenbelt, MD 20771
!                 Phone: 301-614-6161
!                 E-mail: slin@dao.gsfc.nasa.gov
!
! The algorithm is based on the following papers:
 
! 1. Lin, S.-J., and R. B. Rood, 1996: Multidimensional flux form semi-
!    Lagrangian transport schemes. Mon. Wea. Rev., 124, 2046-2070.
!
! 2. Lin, S.-J., W. C. Chao, Y. C. Sud, and G. K. Walker, 1994: A class of
!    the van Leer-type transport schemes and its applications to the moist-
!    ure transport in a General Circulation Model. Mon. Wea. Rev., 122,
!    1575-1593.
!
! 3. Lin, S.-J., and R. B. Rood, 1997: Multidimensional flux form semi-
!    Lagrangian transport schemes- MFCT option. To be submitted.
 
! ======
! INPUT:
! ======
 
! IGD: (horizontal) grid type on which winds are defined.
! IGD = 0  A-Grid  [all variables defined at the same point from south
!                   pole (j=1) to north pole (j=JM) ]

! IGD = 1  GEOS-GCM C-Grid (Max Suarez's center difference dynamical core)

!                                  [North]

!                                   V(i,j)
!                                      |
!                                      |
!                                      |
!                 [WEST]  U(i-1,j)---Q(i,j)---U(i,j) [EAST]
!                                      |
!                                      |
!                                      |
!                                   V(i,j-1)
 
!                                   [South]

!         U(i,   1) is defined at South Pole.
!         V(i,   1) is half grid north of the South Pole.
!         V(i,JM-1) is half grid south of the North Pole.
!
!         V must be defined at j=1 and j=JM-1 if IGD=1
!         V at JM need not be defined.
 
! Q(IM,JM,NL,NC): mixing ratios at current time (t)
! NC: total # of constituents
! IM: first (E-W) dimension; # of Grid intervals in E-W is IM
! JM: 2nd (N-S) dimension;   # of Grid intervals in N-S is JM-1
! NL: 3rd dimension (# of layers); vertical index increases from 1 at
!       the model top to NL near the surface (see fig. below).
!       It is assumed that NL > 5.
!
! PS1(IM,JM): surface pressure at current time (t)
! PS2(IM,JM): surface pressure at mid-time-level (t+NDT/2)
! PS2 is replaced by the predicted PS (at t+NDT) on output.
! Note: surface pressure can have any unit or can be multiplied by any
!       const.
!
! The hybrid ETA-coordinate:
!
! pressure at layer edges are defined as follows:
!
!        p(i,j,k) = AP(k)*PT  +  BP(k)*PS(i,j)          (1)
!
! Where PT is a constant having the same unit as PS.
! AP and BP are unitless constants given at layer edges.
! In all cases  BP(1) = 0., BP(NL+1) = 1.
! The pressure at the model top is PTOP = AP(1)*PT
!
! *********************
! For pure sigma system
! *********************
! AP(k) = 1 for all k, PT = PTOP,
! BP(k) = sige(k) (sigma at edges), PS = Psfc - PTOP, where Psfc
! is the true surface pressure.
!
!                  /////////////////////////////////
!              / \ ------ Model top P=PTOP ---------  AP(1), BP(1)
!               |
!    delp(1)    |  ........... Q(i,j,1) ............
!               |
!     W(k=1)   \ / ---------------------------------  AP(2), BP(2)
!
!
!
!     W(k-1)   / \ ---------------------------------  AP(k), BP(k)
!               |
!    delp(K)    |  ........... Q(i,j,k) ............
!               |
!      W(k)    \ / ---------------------------------  AP(k+1), BP(k+1)
!
!
!
!              / \ ---------------------------------  AP(NL), BP(NL)
!               |
!    delp(NL)   |  ........... Q(i,j,NL) .........
!               |
!     W(NL)=0  \ / -----Earth's surface P=Psfc ------ AP(NL+1), BP(NL+1)
!                 //////////////////////////////////
 
! U(IM,JM,NL) & V(IM,JM,NL):winds (m/s) at mid-time-level (t+NDT/2)
! Note that on return U and V are destroyed.
 
! NDT (integer): time step in seconds (need not be constant during the course of
!      the integration). Suggested value: 30 min. for 4x5, 15 min. for 2x2.5
!      (Lat-Lon) resolution. Smaller values maybe needed if the model
!      has a well-resolved stratosphere and Max(V) > 225 m/s
!
! J1 determines the size of the polar cap:
!    South polar cap edge is located at -90 + (j1-1.5)*180/(JM-1) deg.
!    North polar cap edge is located at  90 - (j1-1.5)*180/(JM-1) deg.
! There are currently only two choices (j1=2 or 3).
! IM must be an even integer if j1 = 2. Recommended value: J1=3.
!
! IORD, JORD, and KORD are integers controlling various options in E-W, N-S,
! and vertical transport, respectively. 
!
!
!  _ORD=
!        1: 1st order upstream scheme (too diffusive, not a real option; it
!           can be used for debugging purposes; this is THE only known "linear"
!           monotonic advection scheme.).
!        2: 2nd order van Leer (full monotonicity constraint;
!           see Lin et al 1994, MWR)
!        3: monotonic PPM* (Collela & Woodward 1984)
!        4: semi-monotonic PPM (same as 3, but overshoots are allowed)
!        5: positive-definite PPM (constraint on the subgrid distribution is
!           only strong enough to prevent generation of negative values;
!           both overshoots & undershootes are possible).
!        6: un-constrained PPM (nearly diffusion free; faster but
!           positivity of the subgrid distribution is not quaranteed. Use
!           this option only when the fields and winds are very smooth or
!           when MFCT=.true.)
!        7: Huynh/Van Leer/Lin full monotonicity constraint
! Only KORD can be set to 7 to enable the use of Huynh's 2nd monotonicity
! constraint for piece-wise parabolic distribution.
!
! *PPM: Piece-wise Parabolic Method
!
! Recommended values:
! IORD=JORD=3 for high horizontal resolution.
! KORD=6 or 7  if MFCT=.true.
! KORD=3 or 7  if MFCT=.false.
!
! The implicit numerical diffusion decreases as _ORD increases.
! DO not use option 4 or 5 for non-positive definite scalars
! (such as Ertel Potential Vorticity).
!
! If numerical diffusion is a problem (particularly at low horizontal
! resolution) then the following setup is recommended:
! IORD=JORD=KORD=6 and MFCT=.true.
!
! AE: Radius of the sphere (meters).
!     Recommended value for the planet earth: 6.371E6
!
! FILL (logical):   flag to do filling for negatives (see note below).
! MFCT (logical):   flag to do a Zalesak-type Multidimensional Flux
!                   correction. It shouldn't be necessary to call the
!                   filling routine when MFCT is true.
!
! Umax: Estimate (upper limit) of the maximum U-wind speed (m/s).
! (225 m/s is a good value for troposphere model; 300 m/s otherwise)
!
! ======
! Output
! ======
!
! Q: the updated mixing ratios at t+NDT (original values are over-written)
! W(;;NL): large-scale vertical mass flux as diagnosed from the hydrostatic
!          relationship. W will have the same unit as PS1 and PS2 (eg, mb).
!          W must be divided by NDT to get the correct mass-flux unit.
!          The vertical Courant number C = W/delp_UPWIND, where delp_UPWIND
!          is the pressure thickness in the "upwind" direction. For example,
!          C(k) = W(k)/delp(k)   if W(k) > 0;
!          C(k) = W(k)/delp(k+1) if W(k) < 0.
!              ( W > 0 is downward, ie, toward surface)
! PS2: predicted PS at t+NDT (original values are over-written)
!
! Memory usage:
! This code is optimized for speed. it requres 18 dynamically allocated
! 3D work arrays (IM,JM,NL) regardless of the value of NC.
! Older versions (version 4 or 4.5) use less memory if NC is small.

! =====
! NOTES:
! =====
!
! This forward-in-time upstream-biased transport scheme degenerates to
! the 2nd order center-in-time center-in-space mass continuity eqn.
! if Q = 1 (constant fields will remain constant). This degeneracy ensures
! that the computed vertical velocity to be identical to GEOS-1 GCM
! for on-line transport.
!
! A larger polar cap is used if j1=3 (recommended for C-Grid winds or when
! winds are noisy near poles).
!
! The user needs to change the parameter Jmax or Kmax if the resolution
! is greater than 0.25 deg in N-S or 500 layers in the vertical direction.
! (this TransPort Core is otherwise resolution independent and can be used
! as a library routine).
 
! PPM is 4th order accurate when grid spacing is uniform (x & y); 3rd
! order accurate for non-uniform grid (vertical sigma coord.).
 
! Time step is limitted only by transport in the meridional direction.
! (the FFSL scheme is not implemented in the meridional direction).
 
! Since only 1-D limiters are applied, negative values could
! potentially be generated when large time step is used and when the
! initial fields contain discontinuities.
! This does not necessarily imply the integration is unstable.
! These negatives are typically very small. A filling algorithm is
! activated if the user set "fill" to be true.
! Alternatively, one can use the MFCT option to enforce monotonicity.

      ! References to F90 modules
      USE GIGC_State_Met_Mod, ONLY : MetState
#if   defined( TOMAS )
      USE TOMAS_MOD,          ONLY : IBINS   !(win, 7/17/09)
      USE TOMAS_TPCORE_MOD,   ONLY : IDFLUX, NBIN, AERMASSDR
#endif

      ! Added to pass C-preprocessor switches (bmy, 3/9/01)
#     include "define.h"

! ****6***0*********0*********0*********0*********0*********0**********72
!clb  PARAMETER (Jmax = 721, kmax = 200)
      PARAMETER (Jmax=1200, kmax=200)
! ****6***0*********0*********0*********0*********0*********0**********72
 
! Input-Output arrays
 
      REAL Q(IM,JM,NL,NC),PS1(IM,JM),PS2(IM,JM),W(IM,JM,NL),
     &     U(IM,JM,NL),V(IM,JM,NL),AP(NL+1),BP(NL+1)
      LOGICAL  ZCROSS, FILL, MFCT, deform

      TYPE(MetState), INTENT(IN) :: State_Met   ! Meteorology State object
 
! Local dynamic arrays
 
      REAL CRX(IM,JM,NL),CRY(IM,JM,NL),delp(IM,JM,NL),delp1(IM,JM,NL),
     &     xmass(IM,JM,NL),ymass(IM,JM,NL),delp2(IM,JM,NL),
     &     DG1(IM),DG2(IM,JM),DPI(IM,JM,NL),qlow(IM,JM,NL),
     &     WK(IM,JM,NL),PU(IM,JM,NL),DQ(IM,JM,NL),
     &     fx(IM+1,JM,NL),fy(IM,JM,NL),fz(IM,JM,NL+1),
     &     qz(IM,JM,NL),Qmax(IM,JM,NL),Qmin(IM,JM,NL)

#if   defined( TOMAS )
! For size-resolved aerosol number and mass (win, 7/17/09)
      REAL  FXNK(IM,JM,NL,IBINS), DQNKX(IM,JM,NL,IBINS)
      REAL  FYNK(IM,JM,NL,IBINS), DQNKY(IM,JM,NL,IBINS)
      REAL  FZNK(IM,JM,NL,IBINS), DQNKZ(IM,JM,NL,IBINS)
#endif

! bey, 6/20/00. for mass-flux diagnostic
      REAL fx1_tp(IM,JM,NL), fy1_tp(IM,JM,NL), fz1_tp(IM,JM,NL)
 
      INTEGER JS(NL),JN(NL)
 
! Local static  arrays
 
      REAL DTDX(Jmax), DTDX5(Jmax), acosp(Jmax),cosp(Jmax),
     &     cose(Jmax), DAP(kmax), DBK(kmax)

      DATA NDT0, NSTEP /0, 0/
      DATA ZCROSS /.true./

! Saved internal variables:
      SAVE DTDY, DTDY5, RCAP, JS0, JN0, IML, DTDX,
     &     DTDX5, acosp, COSP, COSE, DAP,DBK

      ! New variables for TPCORE pressure fixer (bdf, bmy, 10/11/01)
      REAL YMASS_PF(IM,JM,NL), XMASS_PF(IM,JM,NL), TEMP(IM,JM,NL)
      LOGICAL PRESSURE_FIX
      PRESSURE_FIX = .TRUE.

!============================================================================
! Cray NOBOUNDS directive will turn off all subscript bounds checking.
! This directive is also legal on SGI compilers (bmy, 4/24/00)
!DIR$ NOBOUNDS
!============================================================================ 
      deform = .false.
      JM1 = JM -1
      IMH = IM/2
      j2 = JM - j1 + 1
 
      NSTEP = NSTEP + 1  
 
!****6***0*********0*********0*********0*********0*********0**********72
! Initialization
!****6***0*********0*********0*********0*********0*********0**********72
 
      ! Moved further down to be done for each tracer (phs, 30/8/07)
!!      ! For mass flux diagnostics (bey, 6/20/00)
!!      fx1_tp(:,:,:) = 0d0
!!      fy1_tp(:,:,:) = 0d0
!!      fz1_tp(:,:,:) = 0d0
!!      
!!      ! Also need to initialize these arrays, so that the flux diagnostics 
!!      ! will be identical for single or multi processor (bmy, 9/29/00)
!!      fx(:,:,:) = 0d0
!!      fy(:,:,:) = 0d0
!!      fz(:,:,:) = 0d0
!!
      ! Need to initialize these arrays in order to avoid 
      ! floating-point exceptions on Alpha (lyj, bmy, 4/19/02)
      YMASS_PF(:,:,:) = 0d0
      XMASS_PF(:,:,:) = 0d0

      if(NSTEP.eq.1) then

      ! Updated output (bmy, 3/13/03)
      WRITE( 6, '(a)' ) REPEAT( '=', 79 ) 
      WRITE( 6, '(a)' ) 'T P C O R E -- FFSL TransPort Core v. 7.1'
      WRITE( 6, '(a)' )
      WRITE( 6, '(a)' ) 'Originally written by S-J Lin'
      WRITE( 6, '(a)' )
      WRITE( 6, '(a)' )
     & 'Modified for GEOS-CHEM by Isabelle Bey, Brendan Field, and'
      WRITE( 6, '(a)' ) 
     & 'Bob Yantosca, with the addition of flux diagnostics and the'
      WRITE( 6, '(a)' ) 'DYN0 pressure fixer from M. Prather'
      WRITE( 6, '(a)' )
      WRITE( 6, '(a)' ) 'Last Modification Date: 8/22/02'
      WRITE( 6, '(a)' )

#if   ( multitask )
      WRITE( 6, '(a)' ) 'TPCORE was compiled for multitasking'
#if   defined( CRAY )
      WRITE( 6, '(a)' ) 'for CRAY'
#elif defined( SGI_MIPS  )
      WRITE( 6, '(a)' ) 'for SGI Origin/Power Challenge machines'
#elif defined( COMPAQ  )
      WRITE( 6, '(a)' ) 'for COMPAQ/HP RISC Alpha machines'
#elif defined( LINUX_PGI )
      WRITE( 6, '(a)' ) 'for Linux environment w/ PGI compiler'
#elif defined( LINUX_IFORT )
      WRITE( 6, '(a)' ) 'for Linux environment w/ Intel IFORT compiler'
#elif defined( SPARC )
      WRITE( 6, '(a)' ) 'for SUN/Sparc machines'
#elif defined( IBM_AIX )
      WRITE( 6, '(a)' ) 'for IBM/AIX machines'
#endif
#endif

      ! Added output on the first time TPCORE is called (bmy, 10/11/01)
      IF ( PRESSURE_FIX ) THEN
         WRITE( 6, '(a)' ) 
         WRITE( 6, '(a)' ) 'TPCORE PRESSURE FIXER is turned ON!'
      ENDIF
         
      if( MFCT ) then
         WRITE( 6, '(a)' )
         WRITE( 6, '(a)' ) 'MFCT option is on!'
      endif

      ! Updated output (bmy, 4/15/02)
      WRITE( 6, '(a)' )
      WRITE( 6, 2 ) IM, JM, NL, j1
 2    FORMAT( 'IM=  ', i3,1x,'JM=  ', i3,1x,'NL=  ',i3,1x,'J1=  ',i3 )
      
      ! Updated output (bmy, 4/15/02)
      WRITE( 6, 3 ) NC, IORD, JORD, KORD, NDT
 3    FORMAT( 'NC=  ',i3,1x,'IORD=',i3,1x,'JORD=',i3,1x,
     &        'KORD=',i3,1x,'NDT= ',i8)

      if(NL.LT.6) then
        write(6,*) 'stop in module tpcore'
        write(6,*) 'NL must be >=6'
        stop
      endif
 
      if(Jmax.lt.JM .or. Kmax.lt.NL) then
        write(6,*) 'stop in module tpcore'
        write(6,*) 'Jmax or Kmax is too small; see documentation'
        stop
      endif

      DO 5 k=1,NL
      DAP(k) = (AP(k+1) - AP(k))*PT
5     DBK(k) =  BP(k+1) - BP(k)
 
      PI = 4. * ATAN(1.)
      DL = 2.*PI / float(IM)
      DP =    PI / float(JM1)
 
      if(IGD.eq.0) then
! Compute analytic cosine at cell edges
            call cosa(cosp,cose,JM,PI,DP)
      else
! Define cosine consistent with GEOS-GCM (using dycore2.0 or later)
            call cosc(cosp,cose,JM,PI,DP)
      endif
 
      do 15 J=2,JM1
15    acosp(j) = 1./cosp(j)
 
! Inverse of the Scaled polar cap area.
 
      agle = (float(j1)-1.5)*DP
      RCAP  = DP / ( float(IM)*(1.-COS(agle)) )
      acosp(1)  = RCAP
      acosp(JM) = RCAP
      ENDIF
 
      if(NDT0 .ne. NDT) then
      DT   = NDT
      NDT0 = NDT

      CR1  = abs(Umax*DT)/(DL*AE)
      MaxDT = DP*AE / abs(Umax) + 0.5
      
      ! Updated output (bmy, 4/15/02)
      WRITE( 6, '(a)' )
      WRITE(6,*)'Largest time step for max(V)=',Umax,' is ',MaxDT

      ! Bug fix: Now stop the run if NDT is too large.  This will make
      ! sure that we don't violate the Courant limit. (bmy, 11/22/02)
      if(MaxDT .lt. abs(NDT)) then
         write(6,*) 'Warning!!! NDT maybe too large!'
         STOP
      endif

      if(CR1.ge.0.95) then
         JS0 = 0
     	   JN0 = 0
         IML = IM-2
         ZTC = 0.
      else
         ZTC = acos(CR1) * (180./PI)
         JS0 = float(JM1)*(90.-ZTC)/180. + 2
         JS0 = max(JS0, J1+1)
         IML = min(6*JS0/(J1-1)+2, 4*IM/5)
         JN0 = JM-JS0+1
      endif

      ! Updated output (bmy, 4/15/02)
      WRITE( 6, '(''ZTC= '', f13.6)') ZTC
      WRITE( 6, 21 ) JS0, JN0, IML
 21   FORMAT( 'JS=  ',i3,1x, 'JN=  ',i3,1x,'IML= ',i3 )
 
      do 22 J=2,JM1
      DTDX(j)  = DT / ( DL*AE*COSP(J) )
      DTDX5(j) = 0.5*DTDX(j)
22    continue
 
      DTDY  = DT /(AE*DP)
      DTDY5 = 0.5*DTDY
 
      ! Updated output (bmy, 4/15/02)
      WRITE( 6, 23 ) J1, J2
 23   FORMAT( 'J1=  ',i3,1x, 'J2=  ',i3 )

      ! Fancy output to stdout (bmy, 3/13/03)
      WRITE( 6, '(a)' ) REPEAT( '=', 79 )
      ENDIF              ! END INITIALIZATION.
 
!****6***0*********0*********0*********0*********0*********0**********72
! Compute Courant number
!****6***0*********0*********0*********0*********0*********0**********72
 
      if(IGD.eq.0) then
 
! Convert winds on A-Grid to Courant # on C-Grid.

#if   defined( multitask  )
#if   defined( CRAY       ) 
!MIC$ do all shared(NL,im,jm1,jm,U,V,dtdx5,dtdy5,CRX,CRY)
!MIC$* private(i,j,k)
#else
!$OMP PARALLEL DO PRIVATE( I, J, K )
#endif
#endif

      do k=1,NL
      do 46 j=2,JM1
      do 46 i=2,IM
46    CRX(i,j,k) = dtdx5(j)*(U(i,j,k)+U(i-1,j,k))
 
! for i=1
      do 48 j=2,JM1
48    CRX(1,j,k) = dtdx5(j)*(U(1,j,k)+U(IM,j,k))
 
      do 49 j=2,JM
      do 49 i=1,IM
49    CRY(i,j,k) = DTDY5*(V(i,j,k)+V(i,j-1,k))
      enddo
      else
! Convert winds on C-grid to Courant #
! Beware of the index shifting!! (GEOS-GCM)

#if   defined( multitask  )
#if   defined( CRAY       )
!MIC$ do all shared(NL,im,jm1,jm,U,V,dtdx,dtdy,CRX,CRY)
!MIC$* private(i,j,k)
#else
!$OMP PARALLEL DO PRIVATE( I, J, K )
#endif
#endif

      DO 65 k=1,NL
      do 50 j=2,JM1
      do 50 i=2,IM
50    CRX(i,j,k) = dtdx(j)*U(i-1,j,k)
 
      do 55 j=2,JM1
55    CRX(1,j,k) = dtdx(j)*U(IM,j,k)
 
      do 60 j=2,JM
      do 60 i=1,IM
60    CRY(i,j,k) = DTDY*V(i,j-1,k)
65    continue
      endif

      !=================================================================
      ! *****  T P C O R E   P R E S S U R E   F I X E R  *****
      !
      ! Run pressure fixer to fix mass conservation problem.  Pressure 
      ! fixer routines PRESS_FIX, DYN0, PFILTR, LOCFLT, and POLFLT
      ! change the mass fluxes so they become consistant with met field 
      ! pressures.  (bdf, bmy, 10/11/01)
      !
      ! NOTE: The pressure fixer is not 100% perfect; tracer mass will
      !       increase on the order of 0.5%/yr.  However, this is much
      !       better than w/o the pressure fixer, where the mass may
      !       increase by as much as 40%/yr.  (bdf, bmy, 10/22/01)
      !=================================================================
      IF ( PRESSURE_FIX ) THEN

         ! Loop over vertical levels -- 
         ! added parallel loop #if statements (bmy, 10/11/01)
#if   defined( multitask  )
#if   defined( CRAY       )
!MIC$ do all shared(NL,IM,JM,JM1,COSE,XMASS_PF,YMASS_PF,DELP2,CRX,CRY)
!MIC$* private(I,J,K,D5)
#else
!$OMP PARALLEL DO PRIVATE( I, J, K, D5 )
#endif
#endif
         DO K = 1, NL

            ! DELP = pressure thickness: 
            ! the pseudo-density in a hydrostatic system.
            DO J = 1, JM
            DO I = 1, IM
               DELP2(I,J,K) = DAP(K) + DBK(K)*PS2(I,J)
            ENDDO
            ENDDO

            ! calculate mass fluxes for pressure fixer.

            ! N-S component
            DO J = J1, J2+1
               D5 = 0.5 * COSE(J)

               DO I = 1, IM
                  YMASS_PF(I,J,K) = 
     &                 CRY(I,J,K) * D5 * (DELP2(I,J,K)+DELP2(I,J-1,K))
               ENDDO
            ENDDO
 
            ! Enlarged polar cap.
            IF(J1.NE.2) THEN    
               DO I=1,IM
                  YMASS_PF(I,1,K) = 0
                  YMASS_PF(I,JM1+1,K) = 0
               ENDDO
            ENDIF
 
            ! E-W component
            DO J = J1, J2
            DO I =  2, IM
               PU(I,J,K) = 0.5 * (DELP2(I,J,K) + DELP2(I-1,J,K))
            ENDDO
            ENDDO
 
            DO J = J1, J2
               PU(1,J,K) = 0.5 * (DELP2(1,J,K) + DELP2(IM,J,K))
            ENDDO
 
            DO J = J1, J2
            DO I =  1, IM
               XMASS_PF(I,J,K) = PU(I,J,K) * CRX(I,J,K)
            ENDDO
            ENDDO

         ENDDO

         !==============================================================
         ! Call PRESS_FIX to apply the pressure fix to the mass fluxes 
         ! XMASS_PF, YMASS_PF.  PRESS_FIX will call routine DYN0, etc.
         !==============================================================
         CALL PRESS_FIX( XMASS_PF, YMASS_PF, NDT, ACOSP, J1, State_Met )

         ! Loop over vertical levels -- 
         ! added parallel loop #if statements (bmy, 10/11/01)
#if   defined( multitask  )
#if   defined( CRAY       )
!MIC$ do all shared(NL,IM,JM,XMASS_PF,PU,YMASS_PF,DELP2,COSE,CRX,CRY)
!MIC$* private(I,J,K,D5)
#else 
!$OMP PARALLEL DO PRIVATE( I, J, K, D5 )
#endif
#endif
         DO K = 1, NL

            ! Recreate the CRX variable with the new values
            ! of XMASS_PF, which has been adjusted by DYN0
            DO J = J1, J2
            DO I =  1, IM
               CRX(I,J,K) = XMASS_PF(I,J,K) / PU(I,J,K)
            ENDDO
            ENDDO

            ! Recreate the CRY variable with the new values
            ! of YMASS_PF, which has been adjusted by DYN0
            DO J = J1, J2+1
               D5 = 0.5 * COSE(J)

               DO I = 1, IM
                  CRY(I,J,K) = YMASS_PF(I,J,K) /
     &                 ( D5 * ( DELP2(I,J,K) + DELP2(I,J-1,K) ) )
               ENDDO
            ENDDO
         ENDDO
      ENDIF                     

      !=================================================================
      ! End of TPCORE PRESSURE FIXER -- continue as usual
      !=================================================================

!****6***0*********0*********0*********0*********0*********0**********72
! Find JN and JS
!****6***0*********0*********0*********0*********0*********0**********72
 
#if   defined( multitask  )
#if   defined( CRAY       )
!MIC$ do all autoscope shared(JS,JN,CRX,CRY,PS2,U,V,DPI,ymass,delp2,PU)
!MIC$* shared(xmass)
!MIC$* private(i,j,k,sum1,sum2,D5)
#else
!$OMP PARALLEL DO PRIVATE( I, J, K, SUM1, SUM2, D5 )
#endif
#endif

      do 1000 k=1,NL
      JS(k) = j1
      JN(k) = j2
 
      do 111 j=JS0,j1+1,-1
      do 111 i=1,IM
      if(abs(CRX(i,j,k)) .GT. 1.) then
            JS(k) = j
            go to 112
      endif
111   continue
112   continue
 
      do 122 j=JN0,j2-1
      do 122 i=1,IM
      if(abs(CRX(i,j,k)) .GT. 1.) then
            JN(k) = j
            go to 133
      endif
122   continue
133   continue
 
!****6***0*********0*********0*********0*********0*********0**********72
! ***** Compute horizontal mass fluxes *****
!****6***0*********0*********0*********0*********0*********0**********72
 
! delp = pressure thickness: the psudo-density in a hydrostatic system.
      do 30 j=1,JM
      do 30 i=1,IM
30    delp2(i,j,k) = DAP(k) + DBK(k)*PS2(i,j)
 
! N-S componenet
 
      do j=j1,j2+1
      D5 = 0.5 * COSE(j)
      do i=1,IM
      ymass(i,j,k) = CRY(i,j,k)*D5*(delp2(i,j,k) + delp2(i,j-1,k))
      enddo
      enddo
 
      DO 75 j=j1,j2
      DO 75 i=1,IM
75    DPI(i,j,k) = (ymass(i,j,k)-ymass(i,j+1,k)) * acosp(j)
 
      if(j1.ne.2) then           ! Enlarged polar cap.
      do 95 i=1,IM
      DPI(i,  2,k) = 0.
95    DPI(i,JM1,k) = 0.
      endif
 
! Poles
      sum1 = ymass(IM,j1  ,k)
      sum2 = ymass(IM,j2+1,k)
      do 98 i=1,IM-1
      sum1 = sum1 + ymass(i,j1  ,k)
98    sum2 = sum2 + ymass(i,j2+1,k)
 
      sum1 = - sum1 * RCAP
      sum2 =   sum2 * RCAP
      do 100 i=1,IM
      DPI(i, 1,k) = sum1
100   DPI(i,JM,k) = sum2
 
! E-W component
      do j=j1,j2
      do i=2,IM
      PU(i,j,k) = 0.5 * (delp2(i,j,k) + delp2(i-1,j,k))
      enddo
      enddo
 
      do j=j1,j2
      PU(1,j,k) = 0.5 * (delp2(1,j,k) + delp2(IM,j,k))
      enddo
 
      DO 110 j=j1,j2
      DO 110 i=1,IM
110   xmass(i,j,k) = PU(i,j,k)*CRX(i,j,k)
 
      DO 120 j=j1,j2
      DO 120 i=1,IM-1
120   DPI(i,j,k) = DPI(i,j,k) + xmass(i,j,k) - xmass(i+1,j,k)
 
      DO 130 j=j1,j2
130   DPI(IM,j,k) = DPI(IM,j,k) + xmass(IM,j,k) - xmass(1,j,k)
 
!****6***0*********0*********0*********0*********0*********0**********72
! Compute Courant number at cell center
!****6***0*********0*********0*********0*********0*********0**********72
 
      DO 135 j=2,JM1
      do 135 i=1,IM-1
      if(CRX(i,j,k)*CRX(i+1,j,k) .gt. 0.) then
         if(CRX(i,j,k) .gt. 0.) then
         U(i,j,k) = CRX(i,j,k)
         else
         U(i,j,k) = CRX(i+1,j,k)
         endif
      else
         U(i,j,k) = 0.
      endif
135   continue
 
      i=IM
      DO 136 j=2,JM1
      if(CRX(i,j,k)*CRX(1,j,k) .gt. 0.) then
         if(CRX(i,j,k) .gt. 0.) then
         U(i,j,k) = CRX(i,j,k)
         else
         U(i,j,k) = CRX(1,j,k)
         endif
      else
         U(i,j,k) = 0.
      endif
136   continue
 
      do 138 j=2,JM1
      do 138 i=1,IM
      if(CRY(i,j,k)*CRY(i,j+1,k) .gt. 0.) then
         if(CRY(i,j,k) .gt. 0.) then
         V(i,j,k) = CRY(i,j,k)
         else
         V(i,j,k) = CRY(i,j+1,k)
         endif
      else
         V(i,j,k) = 0.
      endif
138   continue
 
      do 139 i=1,IMH
      V(i,     1,k) = 0.5*(CRY(i,2,k)-CRY(i+IMH,2,k))
      V(i+IMH, 1,k) = -V(i,1,k)
      V(i,    JM,k) = 0.5*(CRY(i,JM,k)-CRY(i+IMH,JM1,k))
139   V(i+IMH,JM,k) = -V(i,JM,k)
1000  continue
 
!****6***0*********0*********0*********0*********0*********0**********72
! Compute vertical mass flux (same dimensional unit as PS)
!****6***0*********0*********0*********0*********0*********0**********72
 
! compute total column mass CONVERGENCE.

#if   defined( multitask  )
#if   defined( CRAY       )
!MIC$ do all autoscope shared(im,jm,DPI,PS1,PS2,W,DBK)
!MIC$* shared(DPI,PS1,PS2,W,DBK)
!MIC$* private(i,j,k,DG1)
#else
!$OMP PARALLEL DO PRIVATE( I, J, K, DG1 )
#endif
#endif

      do 395 j=1,jm

      do 320 i=1,IM
320   DG1(i) = DPI(i,j,1)
 
      do 330 k=2,NL
      do 330 i=1,IM
      DG1(i)  = DG1(i) + DPI(i,j,k)
330   continue
 
      do 360 i=1,IM
 
! Compute PS2 (PS at n+1) using the hydrostatic assumption.
! Changes (increases) to surface pressure = total column mass convergence
 
      PS2(i,j)  = PS1(i,j) + DG1(i)
 
! compute vertical mass flux from mass conservation principle.
 
      W(i,j,1) = DPI(i,j,1) - DBK(1)*DG1(i)
      W(i,j,NL) = 0.
360   continue
 
      do 370 k=2,NL-1
      do 370 i=1,IM
      W(i,j,k) = W(i,j,k-1) + DPI(i,j,k) - DBK(k)*DG1(i)
370   continue
395   continue

#if   defined( multitask  )
#if   defined( CRAY       )
!MIC$ do all
!MIC$* shared(deform,NL,im,jm,delp,delp1,delp2,DPI,DAP,DBK,PS1,PS2)
!MIC$* private(i,j,k)
#else
!$OMP PARALLEL DO PRIVATE( I, J, K )
#endif
#endif

      DO 390 k=1,NL

      DO 380 j=1,JM
      DO 380 i=1,IM
      delp1(i,j,k) = DAP(k) + DBK(k)*PS1(i,j)
      delp2(i,j,k) = DAP(k) + DBK(k)*PS2(i,j)
380   delp (i,j,k) = delp1(i,j,k) + DPI(i,j,k)
 
! Check deformation of the flow fields
      if(deform) then

      DO 385 j=1,JM
      DO 385 i=1,IM
      if(delp(i,j,k) .le. 0.) then
c        write(6,*) k,'Noisy wind fields -> delp* is negative!'
c        write(6,*) ' *** Smooth the wind fields or reduce NDT'
         stop
      endif
385   continue
      endif
390   continue

!****6***0*********0*********0*********0*********0*********0**********72
! Do transport one tracer at a time.
!****6***0*********0*********0*********0*********0*********0**********72
 
      DO 5000 IC=1,NC

      ! Moved initialization to 0 here (30/8/07, phs)
      ! For mass flux diagnostics (bey, 6/20/00)
      fx1_tp(:,:,:) = 0d0
      fy1_tp(:,:,:) = 0d0
      fz1_tp(:,:,:) = 0d0
               
      ! Also need to initialize these arrays, so that the flux diagnostics 
      ! will be identical for single or multi processor (bmy, 9/29/00)
      fx(:,:,:) = 0d0
      fy(:,:,:) = 0d0
      fz(:,:,:) = 0d0

#if   defined( TOMAS )
      ! For aerosol mass flux modification (win, 7/17/09)
      FXNK(:,:,:,:) = 0d0
      FYNK(:,:,:,:) = 0d0
      FZNK(:,:,:,:) = 0d0
      DQNKX(:,:,:,:) = 0d0
      DQNKY(:,:,:,:) = 0d0
      DQNKZ(:,:,:,:) = 0d0
#endif


#if   defined( multitask  )
#if   defined( CRAY       )
!MIC$ do all autoscope
!MIC$* shared(q,DQ,delp1,U,V,j1,j2,JS,JN,im,jm,IML,IC,IORD,JORD)
!MIC$* shared(CRX,CRY,PU,xmass,ymass,fx,fy,acosp,rcap,qz)
!MIC$* shared(fx1_tp, fy1_tp)
!MIC$* private(i,j,k,jt,wk,DG2)
#else
!$OMP PARALLEL DO PRIVATE( I, J, K, JT, WK, DG2, IC )
#endif
#endif

      do 2500 k=1,NL
 
      if(j1.ne.2) then
      DO 405 I=1,IM
      q(I,  2,k,IC) = q(I, 1,k,IC)
405   q(I,JM1,k,IC) = q(I,JM,k,IC)
      endif
 
! Initialize DQ
 
      DO 420 j=1,JM
      DO 420 i=1,IM
420   DQ(i,j,k) = q(i,j,k,IC)*delp1(i,j,k)

! E-W advective cross term
      call xadv(IM,JM,j1,j2,q(1,1,k,IC),U(1,1,k),JS(k),JN(k),IML,
     &          wk(1,1,1))
      do 430 j=1,JM
      do 430 i=1,IM
430   wk(i,j,1) = q(i,j,k,IC) + 0.5*wk(i,j,1)
 
! N-S advective cross term
      do 66 j=j1,j2
      do 66 i=1,IM
      jt = float(j) - V(i,j,k)
66    wk(i,j,2) = V(i,j,k) * (q(i,jt,k,IC) - q(i,jt+1,k,IC))
 
      do 77 j=j1,j2
      do 77 i=1,IM
77    wk(i,j,2) = q(i,j,k,IC) + 0.5*wk(i,j,2)

!****6***0*********0*********0*********0*********0*********0**********72
! compute flux in  E-W direction
! Return flux contribution from TPCORE in FX1_TP array (bey, 9/28/00)
      call xtp(IM,JM,IML,j1,j2,JN(k),JS(k),PU(1,1,k),DQ(1,1,k),
     &         wk(1,1,2),CRX(1,1,k),fx(1,1,k),xmass(1,1,k),IORD,
     &         fx1_tp(:,:,k))

! compute flux in  N-S direction
! Return flux contribution from TPCORE in FY1_TP array (bey, 9/28/00)
      call ytp(IM,JM,j1,j2,acosp,RCAP,DQ(1,1,k),wk(1,1,1),
     &         CRY(1,1,k),DG2,ymass(1,1,k),WK(1,1,3),wk(1,1,4),
     &         WK(1,1,5),WK(1,1,6),fy(1,1,k),JORD,
     &         fy1_tp(:,:,k))
!****6***0*********0*********0*********0*********0*********0**********72

      if(ZCROSS) then

! qz is the horizontal advection modified value for input to the
! vertical transport operator FZPPM
! Note: DQ contains only first order upwind contribution.

      do 88 j=1,JM
      do 88 i=1,IM
88    qz(i,j,k) = DQ(i,j,k) / delp(i,j,k)

      else

      do 99 j=1,JM
      do 99 i=1,IM
99    qz(i,j,k) = q(i,j,k,IC)

      endif
 
2500  continue     ! k-loop

!****6***0*********0*********0*********0*********0*********0**********72
! Compute fluxes in the vertical direction
! Return flux contribution from FZPPM in FZ1_TP for ND26 (bey, 9/28/00)
      call FZPPM(qz,fz,IM,JM,NL,DQ,W,delp,KORD,fz1_tp)
!****6***0*********0*********0*********0*********0*********0**********72

      if( MFCT ) then
 
! qlow is the low order "monotonic" solution
 
#if   defined( multitask  )
#if   defined( CRAY       )
!MIC$ do all
!MIC$* shared(NL,im,jm,j1,jm1,qlow,DQ,delp2)
!MIC$* private(i,j,k)
#else
!$OMP PARALLEL DO PRIVATE( I, J, K )
#endif
#endif


      DO k=1,NL

      DO 560 j=1,JM
      DO 560 i=1,IM
560   qlow(i,j,k) = DQ(i,j,k) / delp2(i,j,k)
 
      if(j1.ne.2) then
      DO 561 i=1,IM
      qlow(i,  2,k) = qlow(i, 1,k)
      qlow(i,JM1,k) = qlow(i,JM,k)
561   CONTINUE
      endif

      enddo
!$OMP END PARALLEL DO
 
!****6***0*********0*********0*********0*********0*********0**********72
       call FCT3D(Q(1,1,1,IC),qlow,fx,fy,fz,IM,JM,NL,j1,j2,delp2,
     &           DPI,qz,wk,Qmax,Qmin,DG2,U,V,acosp,RCAP)
! Note: Q is destroyed!!!
!****6***0*********0*********0*********0*********0*********0**********72
      ENDIF
 
! Final update

#if   defined( TOMAS )
!With the addition for TOMAS, I don't know exactly how to make the 
!new code parallel, so c/o this part for now (win, 8/11/09)
c$$$#if   defined( multitask  )
c$$$#if   defined( CRAY       )
c$$$CMIC$ do all autoscope
c$$$CMIC$* private(i,j,k,sum1,sum2)
c$$$#else
c$$$!$OMP PARALLEL DO PRIVATE( I, J, K, SUM1, SUM2 )
c$$$#endif
c$$$#endif 

!Prior to addition for TOMAS aerosols advection (win, 7/17/09)
c$$$      do 101 k=1,NL
c$$$
c$$$      do 425 j=j1,j2
c$$$      do 425 i=1,IM
c$$$      DQ(i,j,k) = DQ(i,j,k) +  fx(i,j,k) - fx(i+1,j,k)
c$$$     &                      + (fy(i,j,k) - fy(i,j+1,k))*acosp(j)
c$$$     &                      +  fz(i,j,k) - fz(i,j,k+1)
c$$$
c$$$425   continue
c$$$
c$$$      sum1 = fy(IM,j1  ,k)
c$$$      sum2 = fy(IM,J2+1,k)
c$$$
c$$$      do i=1,IM-1
c$$$         sum1 = sum1 + fy(i,j1  ,k)
c$$$         sum2 = sum2 + fy(i,J2+1,k)
c$$$      enddo
c$$$ 
c$$$      DQ(1, 1,k) = DQ(1, 1,k) - sum1*RCAP + fz(1, 1,k) - fz(1, 1,k+1)
c$$$      DQ(1,JM,k) = DQ(1,JM,k) + sum2*RCAP + fz(1,JM,k) - fz(1,JM,k+1)
c$$$ 
c$$$      do i=2,IM
c$$$      DQ(i, 1,k) = DQ(1, 1,k)
c$$$      DQ(i,JM,k) = DQ(1,JM,k)
c$$$      enddo
c$$$
c$$$101   continue
!-------ending of original for the modified section (win, 7/17/09)

! Modification for TOMAS aerosol advection (win, 7/17/09)
!Split up the addition of fluxes to DQ array into separate directions;
!E-W, N-S, then vertical.  Then add saving of aerosol number fluxes and 
!aerosol number mixing ratio field and calculation of aerosol mass fluxes
!in between each flux addition to the mixing ratio field.


      IF ( IDFLUX( IC ) /= 0 ) THEN
         CALL AERMASSDR(IC,IM,JM,J1,J2,NL,FXNK(:,:,:,NBIN( IC )),
     &           DQNKZ(:,:,:,NBIN( IC )),fx(1:IM,:,:),DQ(:,:,:),1     )
      ENDIF

      DO K = 1, NL
      DO J = 1, JM
         fx(IM+1, J, K) = fx(1, J, K)
      ENDDO
      ENDDO
         
! Update in E-W direction
      do 425 k=1,NL
      do 425 j=j1,j2
      do 425 i=1,IM
      DQ(i,j,k) = DQ(i,j,k) +  fx(i,j,k) - fx(i+1,j,k)
 425  continue

      IF ( IDFLUX( IC ) /= 0 ) THEN
         CALL AERMASSDR(IC,IM,JM,J1,J2,NL,FYNK(:,:,:,NBIN( IC )),
     &           DQNKX(:,:,:,NBIN( IC )),fy(1:IM,:,:),DQ(:,:,:),2     )
      ENDIF

! Update in N-S direction
      do 101 k=1,NL

      do 426 j=j1,j2
      do 426 i=1,IM
      DQ(i,j,k) = DQ(i,j,k) + (fy(i,j,k) - fy(i,j+1,k))*acosp(j)
426   continue

      sum1 = fy(IM,j1  ,k)
      sum2 = fy(IM,J2+1,k)

      do i=1,IM-1
         sum1 = sum1 + fy(i,j1  ,k)
         sum2 = sum2 + fy(i,J2+1,k)
      enddo
 
      DQ(1, 1,k) = DQ(1, 1,k) - sum1*RCAP 
      DQ(1,JM,k) = DQ(1,JM,k) + sum2*RCAP 
 
      do i=2,IM
      DQ(i, 1,k) = DQ(1, 1,k)
      DQ(i,JM,k) = DQ(1,JM,k)
      enddo

      if(j1.ne.2)then
      do I=1,IM
      DQ(I,  2,k) = DQ(I, 1,k)
      DQ(I,JM1,k) = DQ(I,JM,k)
      enddo
      endif

 101  continue

      IF ( IDFLUX( IC ) /= 0 ) THEN
         CALL AERMASSDR(IC,IM,JM,J1,J2,NL,FZNK(:,:,:,NBIN( IC )),
     &           DQNKY(:,:,:,NBIN( IC )),fz(1:IM,:,:),DQ(:,:,:),3     )
      ENDIF

! Update in the vertical direction
      do 102 k=1,NL

      do 427 j=j1,j2
      do 427 i=1,IM
      DQ(i,j,k) = DQ(i,j,k) +  fz(i,j,k) - fz(i,j,k+1)
 427  continue

      DQ(1, 1,k) = DQ(1, 1,k) + fz(1, 1,k) - fz(1, 1,k+1)
      DQ(1,JM,k) = DQ(1,JM,k) + fz(1,JM,k) - fz(1,JM,k+1)
 
      do i=2,IM
      DQ(i, 1,k) = DQ(1, 1,k)
      DQ(i,JM,k) = DQ(1,JM,k)
      enddo

      if(j1.ne.2)then
      do I=1,IM
      DQ(I,  2,k) = DQ(I, 1,k)
      DQ(I,JM1,k) = DQ(I,JM,k)
      enddo
      endif

 102  continue
!----ending section modified for size-resolved aerosols (win, 7/17/09)

#else

# if   defined( multitask  )
#  if   defined( CRAY       )
!MIC$ do all autoscope
!MIC$* private(i,j,k,sum1,sum2)
#  else
!$OMP PARALLEL DO PRIVATE( I, J, K, SUM1, SUM2 )
#  endif
# endif 

      do 101 k=1,NL

      do 425 j=j1,j2
      do 425 i=1,IM
      DQ(i,j,k) = DQ(i,j,k) +  fx(i,j,k) - fx(i+1,j,k)
     &                      + (fy(i,j,k) - fy(i,j+1,k))*acosp(j)
     &                      +  fz(i,j,k) - fz(i,j,k+1)

425   continue

      sum1 = fy(IM,j1  ,k)
      sum2 = fy(IM,J2+1,k)

      do i=1,IM-1
         sum1 = sum1 + fy(i,j1  ,k)
         sum2 = sum2 + fy(i,J2+1,k)
      enddo
 
      DQ(1, 1,k) = DQ(1, 1,k) - sum1*RCAP + fz(1, 1,k) - fz(1, 1,k+1)
      DQ(1,JM,k) = DQ(1,JM,k) + sum2*RCAP + fz(1,JM,k) - fz(1,JM,k+1)
 
      do i=2,IM
      DQ(i, 1,k) = DQ(1, 1,k)
      DQ(i,JM,k) = DQ(1,JM,k)
      enddo

101   continue
#endif

 
      !=================================================================
      ! bey, 6/20/00. for mass-flux diagnostic
      ! NOTE: DIAG_FLUX is not called within a parallel loop, 
      ! so parallelization can be done within the subroutine
      !=================================================================
      CALL DIAG_FLUX( IC, FX, FX1_TP, FY,  FY1_TP, 
     &                    FZ, FZ1_TP, NDT, ACOSP )

!****6***0*********0*********0*********0*********0*********0**********72
      if(FILL) call qckxyz(DQ,DG2,IM,JM,NL,j1,j2,cosp,acosp,IC,NSTEP)
!****6***0*********0*********0*********0*********0*********0**********72
 
#if   defined( multitask  )
#if   defined( CRAY       )
!MIC$ do all
!MIC$* shared(q,IC,NL,j1,im,jm,jm1,DQ,delp2)
!MIC$* private(i,j,k)
#else
!$OMP PARALLEL DO PRIVATE( I, J, K )
#endif
#endif

      DO k=1,NL

      DO 447 j=1,JM
      DO 447 i=1,IM
447   Q(i,j,k,IC) = DQ(i,j,k) / delp2(i,j,k)
 
      if(j1.ne.2) then
      DO 450 I=1,IM
      Q(I,  2,k,IC) = Q(I, 1,k,IC)
      Q(I,JM1,k,IC) = Q(I,JM,k,IC)
450   CONTINUE
      endif

      enddo     

5000  continue
      RETURN
      END SUBROUTINE TPCORE

!------------------------------------------------------------------------------

      subroutine cosa(cosp,cose,JM,PI,DP)
!****6***0*********0*********0*********0*********0*********0**********72
      REAL cosp(*),cose(*),sine(JM)
!============================================================================
! Cray NOBOUNDS directive will turn off all subscript bounds checking.
! This directive is also legal on SGI compilers (bmy, 4/24/00)
!DIR$ NOBOUNDS
!============================================================================
      do 10 j=2,JM
         ph5  = -0.5*PI + (float(j-1)-0.5)*DP
10    sine(j) = SIN(ph5)
 
      do 80 J=2,JM-1
80    cosp(J) = (sine(j+1)-sine(j))/DP
 
      cosp( 1) = 0.
      cosp(JM) = 0.
 
! Define cosine at edges..

      do 90 j=2,JM
90    cose(j) = 0.5 * (cosp(j-1)+cosp(j))
      cose(1) = cose(2)
      return
      end subroutine cosa

!------------------------------------------------------------------------------

      subroutine cosc(cosp,cose,JNP,PI,DP)
!****6***0*********0*********0*********0*********0*********0**********72
      REAL cosp(*),cose(*)
!============================================================================
! Cray NOBOUNDS directive will turn off all subscript bounds checking.
! This directive is also legal on SGI compilers (bmy, 4/24/00)
!DIR$ NOBOUNDS
!============================================================================
      phi = -0.5*PI
      do 55 j=2,JNP-1
      phi  =  phi + DP
55    cosp(j) = cos(phi)
      cosp(  1) = 0.
      cosp(JNP) = 0.
!
      do 66 j=2,JNP
      cose(j) = 0.5*(cosp(j)+cosp(j-1))
66    CONTINUE
!
      do 77 j=2,JNP-1
      cosp(j) = 0.5*(cose(j)+cose(j+1))
77    CONTINUE
      return
      end subroutine cosc

!------------------------------------------------------------------------------

      subroutine FCT3D(P,plow,fx,fy,fz,im,jm,km,j1,j2,delp,adx,ady,
     &                 wk1,Qmax,Qmin,wkx,CRX,CRY,acosp,RCAP)
!****6***0*********0*********0*********0*********0*********0**********72

      ! Added to pass C-preprocessor switches (bmy, 3/9/01)
#     include "define.h"
 
! MFCT Limiter
! plow: low order solution matrix
! P: current solution matrix

      PARAMETER (esl = 1.E-30)
      REAL P(IM,JM,km),CRX(IM,JM,km),CRY(IM,JM,km),plow(IM,JM,km),
     &     Qmax(IM,JM,km),Qmin(IM,JM,km),acosp(*),delp(im,jm,km),
     &     adx(IM,JM,km),ady(IM,JM,km),fx(IM+1,JM,km),
     &     fy(IM,JM,km),fz(im,jm,km+1),wk1(IM,JM,km),
     &     wkx(im,jm),wkn(im,jm)
!============================================================================
! Cray NOBOUNDS directive will turn off all subscript bounds checking.
! This directive is also legal on SGI compilers (bmy, 4/24/00)
!DIR$ NOBOUNDS
!============================================================================
 
      JM1 = JM-1
 
! Find local min/max of the low-order monotone solution
      call hilo3D(P,im,jm,km,j1,j2,adx,ady,Qmax,Qmin,wkx,wkn)
      call hilo3D(plow,im,jm,km,j1,j2,Qmax,Qmin,wk1,P,wkx,wkn)
! P is destroyed!
 
!     GOTO 123
#if   defined( multitask  )
#if   defined( CRAY       )
!MIC$ do all autoscope
!MIC$* shared(im,j1,j2,km,CRX,CRY,adx,ady,Qmax,Qmin)
!MIC$* private(i,j,k,IT,JT,PS1,PS2,PN1,PN2)
#else
!$OMP PARALLEL DO PRIVATE( I, J, K, IT, JT, PS1, PS2, PN1, PN2 )
#endif
#endif

      DO 1000 k=1,km
      do j=j1,j2
      DO i=1,IM
 
      IT = NINT( float(i) - CRX(i,j,k) )
! Wrap around in E-W
      if(IT .lt. 1) then
            IT = IM + IT
      elseif(IT .GT. IM) then
            IT = IT - IM
      endif
 
      JT = NINT( float(j) - CRY(i,j,k) )
      Qmax(i,j,k) = max(Qmax(i,j,k), adx(IT,JT,k))
      Qmin(i,j,k) = min(Qmin(i,j,k), ady(IT,JT,k))
      enddo
      enddo
 
! Poles:
      PS1 = max(Qmax(1, 1,k), adx(1, 1,k))
      PS2 = min(Qmin(1, 1,k), ady(1, 1,k))
 
      PN1 = max(Qmax(1,JM,k), adx(1,JM,k))
      PN2 = min(Qmin(1,JM,k), ady(1,JM,k))
      DO i=1,IM
      Qmax(i, 1,k) = PS1
      Qmin(i, 1,k) = PS2
 
      Qmax(i,JM,k) = PN1
      Qmin(i,JM,k) = PN2
      enddo
1000  continue
 
123   continue
! Flux Limiter
#if   defined( multitask  )
#if   defined( CRAY       )
!MIC$ do all autoscope
!MIC$* shared(adx,ady,fx,fy,fz,plow,Qmax,Qmin,delp)
!MIC$* private(wkx,wkn)
!MIC$* private(i,j,k,ain,aou,bin,bou,cin,cou,btop,bdon)
#else
!$OMP PARALLEL DO PRIVATE( WKX, WKN, I, J, K, AIN, AOU, 
!$OMP+                     BIN, BOU, CIN, COU, BTOP, BDON )
#endif
#endif

      DO 2000 k=1,km
 
      DO j=j1,j2
      DO i=1,IM
      if(fx(i,j,k) .gt. 0.) then
      Ain = fx(i,j,k)
      Aou = 0.
      else
      Ain = 0.
      Aou = -fx(i,j,k)
      endif
 
      if(fx(i+1,j,k) .gt. 0.) then
      Aou = Aou + fx(i+1,j,k)
      else
      Ain = Ain - fx(i+1,j,k)
      endif
 
      if(fy(i,j,k) .gt. 0.) then
      Bin = fy(i,j,k)
      Bou = 0.
      else
      Bin = 0.
      Bou = -fy(i,j,k)
      endif
 
      if(fy(i,j+1,k) .gt. 0.) then
      Bou = Bou + fy(i,j+1,k)
      else
      Bin = Bin - fy(i,j+1,k)
      endif
 
      if(fz(i,j,k) .gt. 0.) then
      Cin = fz(i,j,k)
      Cou = 0.
      else
      Cin = 0.
      Cou = -fz(i,j,k)
      endif
 
      if(fz(i,j,k+1) .gt. 0.) then
      Cou = Cou + fz(i,j,k+1)
      else
      Cin = Cin - fz(i,j,k+1)
      endif
 
!****6***0*********0*********0*********0*********0*********0**********72
      wkx(i,j) = Ain + Bin*acosp(j) + Cin
      wkn(i,j) = Aou + Bou*acosp(j) + Cou
!****6***0*********0*********0*********0*********0*********0**********72
      enddo
      enddo
 
      DO j=j1,j2
      DO i=1,IM
      adx(i,j,k) = delp(i,j,k)*(Qmax(i,j,k)-plow(i,j,k))/(wkx(i,j)+esl)
      ady(i,j,k) = delp(i,j,k)*(plow(i,j,k)-Qmin(i,j,k))/(wkn(i,j)+esl)
      enddo
      enddo
 
! S Pole
      Ain = 0.
      Aou = 0.
      DO i=1,IM
      if(fy(i,j1,k).gt. 0.) then
           Aou = Aou + fy(i,j1,k)
      else
           Ain = Ain + fy(i,j1,k)
      endif
      enddo
      Ain = -Ain * RCAP
      Aou =  Aou * RCAP
 
! add vertical contribution...
 
      i=1
      j=1
      if(fz(i,j,k) .gt. 0.) then
      Cin = fz(i,j,k)
      Cou = 0.
      else
      Cin = 0.
      Cou = -fz(i,j,k)
      endif
 
      if(fz(i,j,k+1) .gt. 0.) then
      Cou = Cou + fz(i,j,k+1)
      else
      Cin = Cin - fz(i,j,k+1)
      endif
 
!****6***0*********0*********0*********0*********0*********0**********72
      btop = delp(1,1,k)*(Qmax(1,1,k)-plow(1,1,k))/(Ain+Cin+esl)
      bdon = delp(1,1,k)*(plow(1,1,k)-Qmin(1,1,k))/(Aou+Cou+esl)
!****6***0*********0*********0*********0*********0*********0**********72
 
      DO i=1,IM
      adx(i,j,k) = btop
      ady(i,j,k) = bdon
      enddo
! N Pole
      J=JM
      Ain = 0.
      Aou = 0.
      DO i=1,IM
      if(fy(i,j2+1,k).gt. 0.) then
           Ain = Ain + fy(i,j2+1,k)
      else
           Aou = Aou + fy(i,j2+1,k)
      endif
      enddo
      Ain =  Ain * RCAP
      Aou = -Aou * RCAP
 
! add vertical contribution...
 
      i=1
      if(fz(i,j,k) .gt. 0.) then
      Cin = fz(i,j,k)
      Cou = 0.
      else
      Cin = 0.
      Cou = -fz(i,j,k)
      endif
 
      if(fz(i,j,k+1) .gt. 0.) then
      Cou = Cou + fz(i,j,k+1)
      else
      Cin = Cin - fz(i,j,k+1)
      endif
 
!****6***0*********0*********0*********0*********0*********0**********72
      btop = delp(1,j,k)*(Qmax(1,j,k)-plow(1,j,k))/(Ain+Cin+esl)
      bdon = delp(1,j,k)*(plow(1,j,k)-Qmin(1,j,k))/(Aou+Cou+esl)
!****6***0*********0*********0*********0*********0*********0**********72
 
      DO i=1,IM
      adx(i,j,k) = btop
      ady(i,j,k) = bdon
      enddo
 
      if(j1 .ne. 2) then
      DO i=1,IM
! SP
      adx(i,2,k) = adx(i,1,k)
      ady(i,2,k) = ady(i,1,k)
! NP
      adx(i,JM1,k) = adx(i,JM,k)
      ady(i,JM1,k) = ady(i,JM,k)
      enddo
      endif
2000  continue
 
#if   defined( multitask  )
#if   defined( CRAY       )
!MIC$ do all autoscope
!MIC$* shared(fz,adx,ady,im,jm,km)
!MIC$* private(i,j,k)
#else
!$OMP PARALLEL DO PRIVATE( I, J, K )
#endif
#endif

      DO 3000 k=1,km
      DO j=j1,j2
      do i=2,IM
      if(fx(i,j,k) .gt. 0.) then
      fx(i,j,k) = min(1.,ady(i-1,j,k),adx(i,j,k))*fx(i,j,k)
      else
      fx(i,j,k) = min(1.,adx(i-1,j,k),ady(i,j,k))*fx(i,j,k)
      endif
      enddo
      enddo
 
! For i=1
      DO j=j1,j2
      if(fx(1,j,k) .gt. 0.) then
      fx(1,j,k) = min(1.,ady(IM,j,k),adx(1,j,k))*fx(1,j,k)
      else
      fx(1,j,k) = min(1.,adx(IM,j,k),ady(1,j,k))*fx(1,j,k)
      endif
      fx(IM+1,j,k) = fx(1,j,k)
      enddo
 
      do j=j1,j2+1
      do i=1,IM
      if(fy(i,j,k) .gt. 0.) then
        fy(i,j,k) = min(1.,ady(i,j-1,k),adx(i,j,k))*fy(i,j,k)
      else
        fy(i,j,k) = min(1.,adx(i,j-1,k),ady(i,j,k))*fy(i,j,k)
      endif
      enddo
      enddo

      if(k .ne. 1) then
      do j=1,jm
      do i=1,im
      if(fz(i,j,k) .gt. 0.) then
        fz(i,j,k) = min(1.,ady(i,j,k-1),adx(i,j,k))*fz(i,j,k)
      else
        fz(i,j,k) = min(1.,adx(i,j,k-1),ady(i,j,k))*fz(i,j,k)
      endif
      enddo
      enddo
      endif

3000  continue
 
      return
      end subroutine fct3d

!------------------------------------------------------------------------------

      subroutine filew(q,qtmp,IMR,JNP,j1,j2,ipx,tiny)
!****6***0*********0*********0*********0*********0*********0**********72
      REAL q(IMR,*),qtmp(JNP,IMR)
!============================================================================
! Cray NOBOUNDS directive will turn off all subscript bounds checking.
! This directive is also legal on SGI compilers (bmy, 4/24/00)
!DIR$ NOBOUNDS
!============================================================================
!
      ipx = 0
! Copy & swap direction for vectorization.
      do 25 i=1,imr
      do 25 j=j1,j2
25    qtmp(j,i) = q(i,j)
!
      do 55 i=2,imr-1
      do 55 j=j1,j2
      if(qtmp(j,i).lt.0.) then
      ipx =  1
c west
      d0 = max(0.,qtmp(j,i-1))
      d1 = min(-qtmp(j,i),d0)
      qtmp(j,i-1) = qtmp(j,i-1) - d1
      qtmp(j,i) = qtmp(j,i) + d1
c east
      d0 = max(0.,qtmp(j,i+1))
      d2 = min(-qtmp(j,i),d0)
      qtmp(j,i+1) = qtmp(j,i+1) - d2
      qtmp(j,i) = qtmp(j,i) + d2 + tiny
      endif
55    continue
c
      i=1
      do 65 j=j1,j2
      if(qtmp(j,i).lt.0.) then
      ipx =  1
c west
      d0 = max(0.,qtmp(j,imr))
      d1 = min(-qtmp(j,i),d0)
      qtmp(j,imr) = qtmp(j,imr) - d1
      qtmp(j,i) = qtmp(j,i) + d1
c east
      d0 = max(0.,qtmp(j,i+1))
      d2 = min(-qtmp(j,i),d0)
      qtmp(j,i+1) = qtmp(j,i+1) - d2
c
      qtmp(j,i) = qtmp(j,i) + d2 + tiny
      endif
65    continue
      i=IMR
      do 75 j=j1,j2
      if(qtmp(j,i).lt.0.) then
      ipx =  1
c west
      d0 = max(0.,qtmp(j,i-1))
      d1 = min(-qtmp(j,i),d0)
      qtmp(j,i-1) = qtmp(j,i-1) - d1
      qtmp(j,i) = qtmp(j,i) + d1
c east
      d0 = max(0.,qtmp(j,1))
      d2 = min(-qtmp(j,i),d0)
      qtmp(j,1) = qtmp(j,1) - d2
c
      qtmp(j,i) = qtmp(j,i) + d2 + tiny
      endif
75    continue
!
      if(ipx.ne.0) then
      do 85 j=j1,j2
      do 85 i=1,imr
85    q(i,j) = qtmp(j,i)
      else

! Pole
      if(q(1,1).lt.0. or. q(1,JNP).lt.0.) ipx = 1
      endif
      return
      end subroutine filew

!------------------------------------------------------------------------------

      subroutine filns(q,IMR,JNP,j1,j2,cosp,acosp,ipy,tiny)
!****6***0*********0*********0*********0*********0*********0**********72
      REAL q(IMR,*),cosp(*),acosp(*)
      LOGICAL first
      DATA first /.true./
      SAVE cap1
!============================================================================
! Cray NOBOUNDS directive will turn off all subscript bounds checking.
! This directive is also legal on SGI compilers (bmy, 4/24/00)
!DIR$ NOBOUNDS
!============================================================================
!
      if(first) then
      DP = 4.*ATAN(1.)/float(JNP-1)
      cap1 = IMR*(1.-COS((j1-1.5)*DP))/DP
      first = .false.
      endif
!
      ipy = 0
      do 55 j=j1+1,j2-1
      DO 55 i=1,IMR
      IF(q(i,j).LT.0.) THEN
      ipy =  1
      dq  = - q(i,j)*cosp(j)
! North
      dn = q(i,j+1)*cosp(j+1)
      d0 = max(0.,dn)
      d1 = min(dq,d0)
      q(i,j+1) = (dn - d1)*acosp(j+1)
      dq = dq - d1
! South
      ds = q(i,j-1)*cosp(j-1)
      d0 = max(0.,ds)
      d2 = min(dq,d0)
      q(i,j-1) = (ds - d2)*acosp(j-1)
      q(i,j) = (d2 - dq)*acosp(j) + tiny
      endif
55    continue
!
      do i=1,imr
      IF(q(i,j1).LT.0.) THEN
      ipy =  1
      dq  = - q(i,j1)*cosp(j1)
! North
      dn = q(i,j1+1)*cosp(j1+1)
      d0 = max(0.,dn)
      d1 = min(dq,d0)
      q(i,j1+1) = (dn - d1)*acosp(j1+1)
      q(i,j1) = (d1 - dq)*acosp(j1) + tiny
      endif
      enddo
!
      j = j2
      do i=1,imr
      IF(q(i,j).LT.0.) THEN
      ipy =  1
      dq  = - q(i,j)*cosp(j)
! South
      ds = q(i,j-1)*cosp(j-1)
      d0 = max(0.,ds)
      d2 = min(dq,d0)
      q(i,j-1) = (ds - d2)*acosp(j-1)
      q(i,j) = (d2 - dq)*acosp(j) + tiny
      endif
      enddo
!
! Check Poles.
      if(q(1,1).lt.0.) then
      dq = q(1,1)*cap1/float(IMR)*acosp(j1)
      do i=1,imr
      q(i,1) = 0.
      q(i,j1) = q(i,j1) + dq
      if(q(i,j1).lt.0.) ipy = 1
      enddo
      endif
!
      if(q(1,JNP).lt.0.) then
      dq = q(1,JNP)*cap1/float(IMR)*acosp(j2)
      do i=1,imr
      q(i,JNP) = 0.
      q(i,j2) = q(i,j2) + dq
      if(q(i,j2).lt.0.) ipy = 1
      enddo
      endif
!
      return
      end subroutine filns

!------------------------------------------------------------------------------

      subroutine fxppm(IMR,IML,UT,P,DC,fx1,fx2,IORD)
!****6***0*********0*********0*********0*********0*********0**********72
      PARAMETER ( R3 = 1./3., R23 = 2./3. )
      REAL UT(*),fx1(*),P(-IML:IMR+IML+1),DC(-IML:IMR+IML+1)
      REAL AR(0:IMR),AL(0:IMR),A6(0:IMR),fx2(*)
!============================================================================
! Cray NOBOUNDS directive will turn off all subscript bounds checking.
! This directive is also legal on SGI compilers (bmy, 4/24/00)
!DIR$ NOBOUNDS
!============================================================================
!     
      LMT = IORD - 3
!
      DO 10 i=1,IMR
10    AL(i) = 0.5*(p(i-1)+p(i)) + (DC(i-1) - DC(i))*R3
!
      do 20 i=1,IMR-1
20    AR(i) = AL(i+1)
      AR(IMR) = AL(1)
!
      do 30 i=1,IMR
30    A6(i) = 3.*(p(i)+p(i)  - (AL(i)+AR(i)))
!
      if(LMT.LE.2) call lmtppm(DC(1),A6(1),AR(1),AL(1),P(1),IMR,LMT)
!
      AL(0) = AL(IMR)
      AR(0) = AR(IMR)
      A6(0) = A6(IMR)
!
! Abs(UT(i)) < 1
      DO i=1,IMR
      IF(UT(i).GT.0.) then
      fx1(i) = P(i-1)
      fx2(i) = AR(i-1) + 0.5*UT(i)*(AL(i-1) - AR(i-1) +
     &                       A6(i-1)*(1.-R23*UT(i)) )
      else
      fx1(i) = P(i)
      fx2(i) = AL(i) - 0.5*UT(i)*(AR(i) - AL(i) +
     &                     A6(i)*(1.+R23*UT(i)))
      endif
      enddo
!
      DO i=1,IMR
      fx2(i) = fx2(i) - fx1(i)
      enddo
      return
      end subroutine fxppm

!------------------------------------------------------------------------------

      subroutine fyppm(C,P,DC,fy1,fy2,IMR,JNP,j1,j2,A6,AR,AL,JORD)
!****6***0*********0*********0*********0*********0*********0**********72
      PARAMETER ( R3 = 1./3., R23 = 2./3. )
      REAL C(IMR,*),fy1(IMR,*),P(IMR,*),DC(IMR,*),fy2(IMR,JNP)
      REAL AR(IMR,JNP),AL(IMR,JNP),A6(IMR,JNP)
!============================================================================
! Cray NOBOUNDS directive will turn off all subscript bounds checking.
! This directive is also legal on SGI compilers (bmy, 4/24/00)
!DIR$ NOBOUNDS
!============================================================================
!
      IMH = IMR / 2
      JMR = JNP - 1
      j11 = j1-1
      IMJM1 = IMR*(J2-J1+2)
      len   = IMR*(J2-J1+3)
      LMT = JORD - 3
!
      DO 10 i=1,IMR*JMR
      AL(i,2) = 0.5*(p(i,1)+p(i,2)) + (DC(i,1) - DC(i,2))*R3
      AR(i,1) = AL(i,2)
10    CONTINUE
!
! Poles:
!
      DO i=1,IMH
      AL(i,1) = AL(i+IMH,2)
      AL(i+IMH,1) = AL(i,2)
!
      AR(i,JNP) = AR(i+IMH,JMR)
      AR(i+IMH,JNP) = AR(i,JMR)
      enddo
!
      do 30 i=1,len
30    A6(i,j11) = 3.*(p(i,j11)+p(i,j11)  - (AL(i,j11)+AR(i,j11)))
!
      if(LMT.le.2) call lmtppm(DC(1,j11),A6(1,j11),AR(1,j11),
     &                         AL(1,j11),P(1,j11),len,LMT)
!
      DO 140 i=1,IMJM1
      IF(C(i,j1).GT.0.) then
      fy1(i,j1) = P(i,j11)
      fy2(i,j1) = AR(i,j11) + 0.5*C(i,j1)*(AL(i,j11) - AR(i,j11) +
     &                         A6(i,j11)*(1.-R23*C(i,j1)) )
      else
      fy1(i,j1) = P(i,j1)
      fy2(i,j1) = AL(i,j1) - 0.5*C(i,j1)*(AR(i,j1) - AL(i,j1) +
     &                        A6(i,j1)*(1.+R23*C(i,j1)))
      endif
140   continue
c
      DO i=1,IMJM1
      fy2(i,j1) = fy2(i,j1) - fy1(i,j1)
      ENDDO
      return
      end subroutine fyppm

!------------------------------------------------------------------------------

      subroutine FZPPM(P,fz,IMR,JNP,NL,DQ,WZ,delp,KORD,fz1_tp)
!****6***0*********0*********0*********0*********0*********0**********72

      ! Added to pass C-preprocessor switches (bmy, 3/9/01)
#     include "define.h" 

      PARAMETER ( R23 = 2./3., R3 = 1./3.)
      REAL WZ(IMR,JNP,NL),P(IMR,JNP,NL),DQ(IMR,JNP,NL),
     &     fz(IMR,JNP,NL+1),delp(IMR,JNP,NL)
! local 2d arrays
      REAL AR(IMR,NL),AL(IMR,NL),A6(IMR,NL),delq(IMR,NL),DC(IMR,NL)

! bey, 6/20/00. for mass-flux diagnostic
      real fz1_tp(IMR,JNP,NL)

      real lac
!     real x, y, z
!     real median
!     median(x,y,z) = min(max(x,y), max(y,z), max(z,x))
!============================================================================
! Cray NOBOUNDS directive will turn off all subscript bounds checking.
! This directive is also legal on SGI compilers (bmy, 4/24/00)
!DIR$ NOBOUNDS
!============================================================================ 
      km = NL
      km1 = NL-1
      LMT = max(KORD - 3, 0)
 
! find global min/max
 
      ! VMAX1D causes bus errors on SGI.  Replace with F90 intrinsic
      ! functions "MAXVAL" and "MINVAL".  These functions produced
      ! identical results as vmax1d in testing.  "MAXVAL" and "MINVAL"
      ! should also execute more efficiently as well. (bmy, 4/24/00)
      Tmax = MAXVAL( P(:,:,1)  )
      Tmin = MINVAL( P(:,:,1)  )
      Bmax = MAXVAL( P(:,:,NL) )
      Bmin = MINVAL( P(:,:,NL) )

#if   defined( multitask  )
#if   defined( CRAY       )
!MIC$ do all autoscope
!MIC$* shared(LMT,Tmax,Tmin,Bmax,Bmin,JNP,IMR)
!MIC$* shared(fz,DQ,WZ,fz1_tp)
!MIC$* private(i,j,k,c1,c2,tmp,qmax,qmin,A1,A2,d1,d2,qm,dp,c3)
!MIC$* private(cmax,cmin,DC,delq,AR,AL,A6,CM,CP, qmp, lac)
#else
!$OMP PARALLEL DO PRIVATE( I, J, K, C1, C2, TMP, QMAX, QMIN, A1, A2,
!$OMP+                     D1, D2, QM, DP, C3, CMAX, CMIN, DC, DELQ,
!$OMP+                     AR, AL, A6, CM, CP, QMP, LAC )
#endif
#endif

      do 4000 j=1,JNP

      do 500 k=2,km
      do 500 i=1,IMR
500   A6(i,k) = delp(i,j,k-1) + delp(i,j,k)

      do 1000 k=1,km1
      do 1000 i=1,IMR
1000  delq(i,k) = P(i,j,k+1) - P(i,j,k)
 
      DO 1220 k=2,km1
      DO 1220 I=1,IMR
      c1 = (delp(i,j,k-1)+0.5*delp(i,j,k))/A6(i,k+1)
      c2 = (delp(i,j,k+1)+0.5*delp(i,j,k))/A6(i,k)
      tmp = delp(i,j,k)*(c1*delq(i,k) + c2*delq(i,k-1))
     &     / (A6(i,k)+delp(i,j,k+1))
      Qmax = max(P(i,j,k-1),P(i,j,k),P(i,j,k+1)) - P(i,j,k)
      Qmin = P(i,j,k) - min(P(i,j,k-1),P(i,j,k),P(i,j,k+1))
      DC(i,k) = sign(min(abs(tmp),Qmax,Qmin), tmp)
1220  CONTINUE
 
!****6***0*********0*********0*********0*********0*********0**********72
! Compute the first guess at cell interface
! First guesses are required to be continuous.
!****6***0*********0*********0*********0*********0*********0**********72
 
! Interior.
 
      DO 12 k=3,km1
      DO 12 i=1,IMR
      c1 = delq(i,k-1)*delp(i,j,k-1) / A6(i,k)
      A1 = A6(i,k-1) / (A6(i,k) + delp(i,j,k-1))
      A2 = A6(i,k+1) / (A6(i,k) + delp(i,j,k))
      AL(i,k) = P(i,j,k-1) + c1 + 2./(A6(i,k-1)+A6(i,k+1)) *
     &          ( delp(i,j,k  )*(c1*(A1 - A2)+A2*DC(i,k-1)) -
     &                          delp(i,j,k-1)*A1*DC(i,k  ) )
12    CONTINUE
 
! Area preserving cubic with 2nd deriv. = 0 at the boundaries
! Top
      DO 10 i=1,IMR
      d1 = delp(i,j,1)
      d2 = delp(i,j,2)
      qm = (d2*P(i,j,1)+d1*P(i,j,2)) / (d1+d2)
      dp = 2.*(P(i,j,2)-P(i,j,1)) / (d1+d2)
      c1 = 4.*(AL(i,3)-qm-d2*dp) / ( d2*(2.*d2*d2+d1*(d2+3.*d1)) )
      c3 = dp - 0.5*c1*(d2*(5.*d1+d2)-3.*d1**2)
      AL(i,2) = qm - 0.25*c1*d1*d2*(d2+3.*d1)
      AL(i,1) = d1*(2.*c1*d1**2-c3) + AL(i,2)
      DC(i,1) =  P(i,j,1) - AL(i,1)
! No over- and undershoot condition
      AL(i,1) = max(Tmin,AL(i,1))
      AL(i,1) = min(Tmax,AL(i,1))
      Cmax = max(P(i,j,1), P(i,j,2))
      Cmin = min(P(i,j,1), P(i,j,2))
      AL(i,2) = max(Cmin,AL(i,2))
      AL(i,2) = min(Cmax,AL(i,2))
10    continue
 
! Bottom
      DO 15 i=1,IMR
      d1 = delp(i,j,km )
      d2 = delp(i,j,km1)
      qm = (d2*P(i,j,km)+d1*P(i,j,km1)) / (d1+d2)
      dp = 2.*(P(i,j,km1)-P(i,j,km)) / (d1+d2)
	c1 = 4.*(AL(i,km1)-qm-d2*dp) / (d2*(2.*d2*d2+d1*(d2+3.*d1)))
	c3 = dp - 0.5*c1*(d2*(5.*d1+d2)-3.*d1**2)
      AL(i,km) = qm - 0.25*c1*d1*d2*(d2+3.*d1)
      AR(i,km) = d1*(2.*c1*d1**2-c3) + AL(i,km)
      DC(i,km) = AR(i,km) -  P(i,j,km)
! No over- and undershoot condition
      Cmax = max(P(i,j,km), P(i,j,km1))
      Cmin = min(P(i,j,km), P(i,j,km1))
      AL(i,km) = max(Cmin,AL(i,km))
      AL(i,km) = min(Cmax,AL(i,km))
      AR(i,km) = max(Bmin,AR(i,km))
      AR(i,km) = min(Bmax,AR(i,km))
15    continue

      do 20 k=1,km1
      do 20 i=1,IMR
      AR(i,k) = AL(i,k+1)
20    continue
 
! f(s) = AL + s*[(AR-AL) + A6*(1-s)]         ( 0 <= s  <= 1 )
! Top 2 layers
      do k=1,2
         do i=1,IMR
         A6(i,k) = 3.*(P(i,j,k)+P(i,j,k) - (AL(i,k)+AR(i,k)))
         enddo
      call lmtppm(DC(1,k),A6(1,k),AR(1,k),AL(1,k),P(1,j,k),
     &            IMR,0)
      enddo

! Interior.
      if(LMT.LE.2) then
      do k=3,NL-2
         do i=1,IMR
         A6(i,k) = 3.*(P(i,j,k)+P(i,j,k) - (AL(i,k)+AR(i,k)))
         enddo
      call lmtppm(DC(1,k),A6(1,k),AR(1,k),AL(1,k),P(1,j,k),
     &            IMR,LMT)
      enddo

      elseif(LMT .eq. 4) then

!****6***0*********0*********0*********0*********0*********0**********72
! Huynh's 2nd constraint
!****6***0*********0*********0*********0*********0*********0**********72

      do k=2, NL-1
         do i=1,imr
            DC(i,k) = delq(i,k) - delq(i,k-1)
         enddo
      enddo

      do  k=3, NL-2
         do  i=1, imr
! Right edges
         qmp   = P(i,j,k)                 + 2.0*delq(i,k-1)
         lac   = P(i,j,k) + 1.5*DC(i,k-1) + 0.5*delq(i,k-1)
         qmin  = min(P(i,j,k), qmp, lac)
         qmax  = max(P(i,j,k), qmp, lac)
!        AR(i,k) = median(AR(i,k), qmin, qmax)
         AR(i,k) = min(max(AR(i,k), qmin), qmax)
! Left  edges
         qmp   = P(i,j,k)                 - 2.0*delq(i,k)
         lac   = P(i,j,k) + 1.5*DC(i,k+1) - 0.5*delq(i,k)
         qmin  = min(P(i,j,k), qmp, lac)
         qmax  = max(P(i,j,k), qmp, lac)
!        AL(i,k) = median(AL(i,k), qmin, qmax)
         AL(i,k) = min(max(AL(i,k), qmin), qmax)
! Recompute A6
         A6(i,k) = 3.*(2.*P(i,j,k) - (AR(i,k)+AL(i,k)))
         enddo
      enddo
      endif

! Bottom 2 layers
      do k=NL-1,NL
         do i=1,IMR
         A6(i,k) = 3.*(P(i,j,k)+P(i,j,k) - (AL(i,k)+AR(i,k)))
         enddo
      call lmtppm(DC(1,k),A6(1,k),AR(1,k),AL(1,k),P(1,j,k),
     &            IMR,0)
      enddo
 
      DO 140 k=2,NL
      DO 140 i=1,IMR
      IF(WZ(i,j,k-1).GT.0.) then
             CM = WZ(i,j,k-1) / delp(i,j,k-1)
      DC(i,k) = P(i,j,k-1)
      fz(i,j,k) = AR(i,k-1)+0.5*CM*(AL(i,k-1)-AR(i,k-1)+
     &                      A6(i,k-1)*(1.-R23*CM))
      else
             CP = WZ(i,j,k-1) / delp(i,j,k)
      DC(i,k) = P(i,j,k)
      fz(i,j,k) = AL(i,k)+0.5*CP*(AL(i,k)-AR(i,k)-
     &                    A6(i,k)*(1.+R23*CP))
      endif
140   continue
 
      DO 250 k=2,NL
      DO 250 i=1,IMR
      fz(i,j,k) = WZ(i,j,k-1) * (fz(i,j,k) - DC(i,k))
      DC(i,k) = WZ(i,j,k-1) * DC(i,k)
250   continue

      do 350 i=1,IMR
      fz(i,j,   1) = 0.
      fz(i,j,NL+1) = 0.
      DQ(i,j, 1) = DQ(i,j, 1) - DC(i, 2)
      DQ(i,j,NL) = DQ(i,j,NL) + DC(i,NL)
      fz1_tp(i,j,1)  = 0.        ! PHS
      fz1_tp(i,j,NL) = DC(i,NL)  ! PHS - flux b/w 1st and second layer
350   continue
 
!-----------------------------------------------------------------------------
! bey, 6/20/00. for mass-flux diagnostic, loop had to be extended 
!      do 360 k=2,km1
!      do 360 i=1,IMR
!360   DQ(i,j,k) = DQ(i,j,k) + DC(i,k) - DC(i,k+1)
!-----------------------------------------------------------------------------
      do k=2,km1
      do i=1,IMR
         DQ(i,j,k) = DQ(i,j,k) + DC(i,k) - DC(i,k+1)

         ! bey, 6/20/00. for mass-flux diagnostic
         fz1_tp(i,j,k) = DC(i,k)
      enddo
      enddo
 
#if   defined( TOMAS )
      ! Extend the loop to k=NL (win, 7/17/09)
      do i=1,IMR
         fz1_tp(i,j,NL) = DC(i,NL)
      enddo
#endif
 
4000  continue
      return
      end subroutine fzppm


!------------------------------------------------------------------------------

      subroutine hilo(q,im,jm,j1,j2,qmax,qmin,bt,bd)
!****6***0*********0*********0*********0*********0*********0**********72
      REAL q(IM,JM),Qmax(IM,JM),Qmin(IM,JM),bt(IM,*),bd(IM,*)
!============================================================================
! Cray NOBOUNDS directive will turn off all subscript bounds checking.
! This directive is also legal on SGI compilers (bmy, 4/24/00)
!DIR$ NOBOUNDS
!============================================================================ 
! y-sweep
      DO j=j1,j2
      DO i=1,IM
      bt(i,j) = max(q(i,j-1),q(i,j),q(i,j+1))
      bd(i,j) = min(q(i,j-1),q(i,j),q(i,j+1))
      enddo
      enddo
!
! x-sweep
      IM1 = IM-1
      DO j=j1,j2
      DO i=2,IM1
      Qmax(i,j) = max(bt(i-1,j),bt(i,j),bt(i+1,j))
      Qmin(i,j) = min(bd(i-1,j),bd(i,j),bd(i+1,j))
      enddo
      enddo
!
      DO j=j1,j2
!     i = 1
      Qmax(1,j) = max(bt(IM,j),bt(1,j),bt(2,j))
      Qmin(1,j) = min(bd(IM,j),bd(1,j),bd(2,j))
!     i = IM
      Qmax(IM,j) = max(bt(IM1,j),bt(IM,j),bt(1,j))
      Qmin(IM,j) = min(bd(IM1,j),bd(IM,j),bd(1,j))
      enddo
!
! N. Pole:
      Pmax = q(1,JM)
      Pmin = q(1,JM)
      do i=1,IM
      if(q(i,j2) .gt. Pmax) then
            Pmax = q(i,j2)
      elseif(q(i,j2) .lt. Pmin) then
            Pmin = q(i,j2)
      endif
      enddo
!
      do i=1,IM
      Qmax(i,JM) = Pmax
      Qmin(i,JM) = Pmin
      enddo
!
! S. Pole:
      Pmax = q(1,1)
      Pmin = q(1,1)
      do i=1,IM
      if(q(i,j1) .gt. Pmax) then
            Pmax = q(i,j1)
      elseif(q(i,j1) .lt. Pmin) then
            Pmin = q(i,j1)
      endif
      enddo
!
      do i=1,IM
      Qmax(i,1) = Pmax
      Qmin(i,1) = Pmin
      enddo
!
      if(j1 .ne. 2) then
      JM1 = JM-1
      do i=1,IM
      Qmax(i,2) = Qmax(i,1)
      Qmin(i,2) = Qmin(i,1)
!
      Qmax(i,JM1) = Qmax(i,JM)
      Qmin(i,JM1) = Qmin(i,JM)
      enddo
      endif
      return
      end subroutine hilo

!------------------------------------------------------------------------------

      subroutine hilo3D(P,im,jm,km,j1,j2,Pmax,Pmin,Qmax,Qmin,bt,bd)
!****6***0*********0*********0*********0*********0*********0**********72

      ! Added to pass C-preprocessor switches (bmy, 3/9/01)
#     include "define.h" 

      REAL P(IM,JM,km),Pmax(IM,JM,km),Pmin(IM,JM,km),
     &     Qmax(IM,JM,km),Qmin(IM,JM,km),bt(im,jm),bd(im,jm)
!============================================================================
! Cray NOBOUNDS directive will turn off all subscript bounds checking.
! This directive is also legal on SGI compilers (bmy, 4/24/00)
!DIR$ NOBOUNDS
!============================================================================ 
#if   defined( multitask  )
#if   defined( CRAY       )
!MIC$ do all autoscope
!MIC$* shared(P,Qmax,Qmin,im,jm,j1,j2)
!MIC$* private(k,bt,bd)
#else
!$OMP PARALLEL DO PRIVATE( K, BT, BD )
#endif
#endif

      DO 1000 k=1,km
      call hilo(P(1,1,k),im,jm,j1,j2,Qmax(1,1,k),Qmin(1,1,k),bt,bd)
1000  continue
 
      km1 = km-1
      km2 = km-2

#if   defined( multitask  )
#if   defined( CRAY       )
!MIC$ do all autoscope
!MIC$* shared(Pmax,Pmin,Qmax,Qmin,im,jm,km,km1,km2)
!MIC$* private(i,j)
#else
!$OMP PARALLEL DO PRIVATE( I, J )
#endif
#endif
 
      DO 2000 j=1,jm
      DO 2000 i=1,im
! k=1 and k=km
      Pmax(i,j, 1) = max(Qmax(i,j,  2),Qmax(i,j, 1))
      Pmin(i,j, 1) = min(Qmin(i,j,  2),Qmin(i,j, 1))
      Pmax(i,j,km) = max(Qmax(i,j,km1),Qmax(i,j,km))
      Pmin(i,j,km) = min(Qmin(i,j,km1),Qmin(i,j,km))
! k=2 and k=km1
      Pmax(i,j,  2) = max(Qmax(i,j,  3),Pmax(i,j, 1))
      Pmin(i,j,  2) = min(Qmin(i,j,  3),Pmin(i,j, 1))
      Pmax(i,j,km1) = max(Qmax(i,j,km2),Pmax(i,j,km))
      Pmin(i,j,km1) = min(Qmin(i,j,km2),Pmin(i,j,km))
2000  continue

#if   defined( multitask  )
#if   defined( CRAY       )
!MIC$ do all autoscope
!MIC$* shared(Pmax,Pmin,Qmax,Qmin,im,jm,km,km1,km2)
!MIC$* private(i,j,k)
#else
!$OMP PARALLEL DO PRIVATE( I, J, K )
#endif
#endif

      DO 3000 k=3,km2
      DO 3000 j=1,jm
      DO 3000 i=1,im
      Pmax(i,j,k) = max(Qmax(i,j,k-1),Qmax(i,j,k),Qmax(i,j,k+1))
      Pmin(i,j,k) = min(Qmin(i,j,k-1),Qmin(i,j,k),Qmin(i,j,k+1))
3000  continue
      return
      end subroutine hilo3D

!------------------------------------------------------------------------------

      subroutine lmtppm(DC,A6,AR,AL,P,IM,LMT)
!****6***0*********0*********0*********0*********0*********0**********72
!
! A6 =  CURVATURE OF THE TEST PARABOLA
! AR =  RIGHT EDGE VALUE OF THE TEST PARABOLA
! AL =  LEFT  EDGE VALUE OF THE TEST PARABOLA
! DC =  0.5 * MISMATCH
! P  =  CELL-AVERAGED VALUE
! IM =  VECTOR LENGTH
!
! OPTIONS:
!
! LMT = 0: FULL MONOTONICITY
! LMT = 1: SEMI-MONOTONIC CONSTRAINT (NO UNDERSHOOTS)
! LMT = 2: POSITIVE-DEFINITE CONSTRAINT
!
      PARAMETER ( R12 = 1./12. )
      REAL A6(IM),AR(IM),AL(IM),P(IM),DC(IM)
!============================================================================
! Cray NOBOUNDS directive will turn off all subscript bounds checking.
! This directive is also legal on SGI compilers (bmy, 4/24/00)
!DIR$ NOBOUNDS
!============================================================================
!
      if(LMT.eq.0) then
! Full constraint
      do 100 i=1,IM
      if(DC(i).eq.0.) then
            AR(i) = p(i)
            AL(i) = p(i)
            A6(i) = 0.
      else
      da1  = AR(i) - AL(i)
      da2  = da1**2
      A6DA = A6(i)*da1
      if(A6DA .lt. -da2) then
            A6(i) = 3.*(AL(i)-p(i))
            AR(i) = AL(i) - A6(i)
      elseif(A6DA .gt. da2) then
            A6(i) = 3.*(AR(i)-p(i))
            AL(i) = AR(i) - A6(i)
      endif
      endif
100   continue
      elseif(LMT.eq.1) then
! Semi-monotonic constraint
      do 150 i=1,IM
      if(abs(AR(i)-AL(i)) .GE. -A6(i)) go to 150
      if(p(i).lt.AR(i) .and. p(i).lt.AL(i)) then
            AR(i) = p(i)
            AL(i) = p(i)
            A6(i) = 0.
      elseif(AR(i) .gt. AL(i)) then
            A6(i) = 3.*(AL(i)-p(i))
            AR(i) = AL(i) - A6(i)
      else
            A6(i) = 3.*(AR(i)-p(i))
            AL(i) = AR(i) - A6(i)
      endif
150   continue
      elseif(LMT.eq.2) then
      do 250 i=1,IM
      if(abs(AR(i)-AL(i)) .GE. -A6(i)) go to 250
      fmin = p(i) + 0.25*(AR(i)-AL(i))**2/A6(i) + A6(i)*R12
      if(fmin.ge.0.) go to 250
      if(p(i).lt.AR(i) .and. p(i).lt.AL(i)) then
            AR(i) = p(i)
            AL(i) = p(i)
            A6(i) = 0.
      elseif(AR(i) .gt. AL(i)) then
            A6(i) = 3.*(AL(i)-p(i))
            AR(i) = AL(i) - A6(i)
      else
            A6(i) = 3.*(AR(i)-p(i))
            AL(i) = AR(i) - A6(i)
      endif
250   continue
      endif
      return
      end subroutine lmtppm

!------------------------------------------------------------------------------

      SUBROUTINE qckxyz(Q,qtmp,IMR,JNP,NLAY,j1,j2,cosp,acosp,IC,NSTEP)
!****6***0*********0*********0*********0*********0*********0**********72

      ! Added to pass C-preprocessor switches (bmy, 3/9/01)
#     include "define.h" 
 
      PARAMETER ( tiny = 1.E-30 )
      PARAMETER ( kmax = 200 )
      REAL Q(IMR,JNP,NLAY),qtmp(IMR,JNP),cosp(*),acosp(*)
      integer IP(kmax)
!============================================================================
! Cray NOBOUNDS directive will turn off all subscript bounds checking.
! This directive is also legal on SGI compilers (bmy, 4/24/00)
!DIR$ NOBOUNDS
!============================================================================ 
      NLM1 = NLAY-1

! Do horizontal filling.

#if   defined( multitask  )
#if   defined( CRAY       )
!MIC$ do all autoscope
!MIC$* private(i,j,L,qtmp)
#else
!$OMP PARALLEL DO PRIVATE( I, J, L, QTMP )
#endif
#endif

      do 1000 L=1,NLAY
         call filns(q(1,1,L),IMR,JNP,j1,j2,cosp,acosp,ip(L),tiny)
         if(ip(L).ne.0) 
     &   call filew(q(1,1,L),qtmp,IMR,JNP,j1,j2,ip(L),tiny)
1000  continue

      ipz = 0
      do L=1,NLAY
      if(ip(L) .ne. 0) then
         ipz = L
         go to 111
      endif
      enddo
      return

111   continue

      if(ipz .eq. 0) return

      if(ipz .eq. 1) then
         lpz = 2
      else
         lpz = ipz
      endif

! Do vertical filling.

#if   defined( multitask  )
#if   defined( CRAY       )
!MIC$ do all autoscope
!MIC$* private(i,j,L,qup,qly,dup)
#else 
!$OMP PARALLEL DO PRIVATE( I, J, L, QUP, QLY, DUP )
#endif
#endif

      do 2000 j=j1,j2

      if(ipz .eq. 1) then
! Top layer
      do i=1,IMR
      if(Q(i,j,1).LT.0.) then
          Q(i,j,2) = Q(i,j,2) + Q(i,j,1)
          Q(i,j,1) = 0.
      endif
      enddo
      endif
 
      DO 225 L = lpz,NLM1
      do i=1,IMR
      IF( Q(i,j,L).LT.0.) THEN
! From above
          qup =  Q(i,j,L-1)
          qly = -Q(i,j,L)
          dup  = min(qly,qup)
          Q(i,j,L-1) = qup - dup
          Q(i,j,L  ) = dup-qly
! Below
          Q(i,j,L+1) = Q(i,j,L+1) + Q(i,j,L)
          Q(i,j,L)   = 0.
      ENDIF
      enddo
225   CONTINUE
 
! BOTTOM LAYER
      L = NLAY
      do i=1,IMR
      IF( Q(i,j,L).LT.0.) THEN
 
! From above
 
          qup = Q(i,j,NLM1)
          qly = -Q(i,j,L)
          dup = min(qly,qup)
          Q(i,j,NLM1) = qup - dup

! From "below" the surface.
          Q(i,j,L) = 0.
      ENDIF
      enddo
2000  continue

      RETURN
      END SUBROUTINE qckxyz

!------------------------------------------------------------------------------

      subroutine xadv(IMR,JNP,j1,j2,p,UA,JS,JN,IML,adx)
!****6***0*********0*********0*********0*********0*********0**********72
      REAL p(IMR,JNP),adx(IMR,JNP),qtmp(-IMR:IMR+IMR),UA(IMR,JNP)
!============================================================================
! Cray NOBOUNDS directive will turn off all subscript bounds checking.
! This directive is also legal on SGI compilers (bmy, 4/24/00)
!DIR$ NOBOUNDS
!============================================================================
!     
      JMR = JNP-1
      do 1309 j=j1,j2
      if(J.GT.JS  .and. J.LT.JN) GO TO 1309
!
      do i=1,IMR
      qtmp(i) = p(i,j)
      enddo
!
      do i=-IML,0
      qtmp(i)       = p(IMR+i,j)
      qtmp(IMR+1-i) = p(1-i,j)
      enddo
!
      DO i=1,IMR
      iu = UA(i,j)
      ru = UA(i,j) - iu
      iiu = i-iu
      if(UA(i,j).GE.0.) then
      adx(i,j) = qtmp(iiu)+ru*(qtmp(iiu-1)-qtmp(iiu))
      else
      adx(i,j) = qtmp(iiu)+ru*(qtmp(iiu)-qtmp(iiu+1))
      endif
      enddo
 
      do i=1,IMR
      adx(i,j) = adx(i,j) - p(i,j)
      enddo
1309  continue
 
! Eulerian upwind
 
      do j=JS+1,JN-1
!
      do i=1,IMR
      qtmp(i) = p(i,j)
      enddo
!
      qtmp(0)     = p(IMR,J)
      qtmp(IMR+1) = p(1,J)
!
      DO i=1,IMR
      IP = i - UA(i,j)
      adx(i,j) = UA(i,j)*(qtmp(ip)-qtmp(ip+1))
      enddo
      enddo
!
      if(j1.ne.2) then
      do i=1,IMR
      adx(i,  2) = 0.
      adx(i,JMR) = 0.
      enddo
      endif

! set cross term due to x-adv at the poles to zero.
      do i=1,IMR
      adx(i,  1) = 0.
      adx(i,JNP) = 0.
      enddo
      return
      end subroutine xadv

!------------------------------------------------------------------------------

      subroutine xmist(IMR,IML,P,DC)
!****6***0*********0*********0*********0*********0*********0**********72
      REAL P(-IML:IMR+1+IML),DC(-IML:IMR+1+IML)
!============================================================================
! Cray NOBOUNDS directive will turn off all subscript bounds checking.
! This directive is also legal on SGI compilers (bmy, 4/24/00)
!DIR$ NOBOUNDS
!============================================================================
!
! 2nd order version.
!
      do 10  i=1,IMR
      tmp = 0.25*(p(i+1) - p(i-1))
      Pmax = max(P(i-1), p(i), p(i+1)) - p(i)
      Pmin = p(i) - min(P(i-1), p(i), p(i+1))
10    DC(i) = sign(min(abs(tmp),Pmax,Pmin), tmp)
      return
      end subroutine xmist

!------------------------------------------------------------------------------

      subroutine xtp(im,jm,IML,j1,j2,JN,JS,PU,DQ,q,C,fx2,xmass,IORD,
     &  fx1_tp)
!****6***0*********0*********0*********0*********0*********0**********72

      REAL C(im,*),DC(-IML:im+IML+1),xmass(im,jm),
     &     fx1(im+1),DQ(im,jm),qtmp(-IML:im+1+IML)
      REAL PU(im,jm),q(im,jm)
      real fx2(im+1,jm)
      INTEGER isave(im)
! bey, 6/20/00. for mass-flux diagnostic
      real fx1_tp(im,jm)

!============================================================================
! Cray NOBOUNDS directive will turn off all subscript bounds checking.
! This directive is also legal on SGI compilers (bmy, 4/24/00)
!DIR$ NOBOUNDS
!============================================================================
!
      IMP = im + 1
!
! van Leer at high latitudes
      jvan = max(1,jm/20)
      j1vl = j1+jvan
      j2vl = j2-jvan
!
      do 1310 j=j1,j2
!
      do i=1,im
         qtmp(i) = q(i,j)
      enddo
!
      if(j.ge.JN .or. j.le.JS) goto 2222
!****6***0*********0*********0*********0*********0*********0**********72
! *** Eulerian ***
!****6***0*********0*********0*********0*********0*********0**********72
!
      qtmp(0)     = q(im,J)
      qtmp(-1)    = q(im-1,J)
      qtmp(IMP)   = q(1,J)
      qtmp(IMP+1) = q(2,J)
 
      IF(IORD.eq.1 .or. j.eq.j1. or. j.eq.j2) THEN
      do i=1,im
         iu = float(i) - c(i,j)
         fx1(i) = qtmp(iu)
      enddo
 
! Zero high order contribution
      DO i=1,im
         fx2(i,j) = 0.
      enddo
      ELSE
      call xmist(im,IML,Qtmp,DC)
      DC(0) = DC(im)
!
      if(IORD.eq.2 .or. j.le.j1vl .or. j.ge.j2vl) then
      DO i=1,im
         iu = float(i) - c(i,j)
         fx1(i  ) = qtmp(iu)
         fx2(i,j) = DC(iu)*(sign(1.,c(i,j))-c(i,j))
      enddo
      else
      call fxppm(im,IML,C(1,j),Qtmp,DC,fx1,fx2(1,j),IORD)
      endif
!
      ENDIF
!
      DO i=1,im
         fx1(i  ) = fx1(i  )*xmass(i,j)
         fx2(i,j) = fx2(i,j)*xmass(i,j)
      enddo
!
      goto 1309
!
!****6***0*********0*********0*********0*********0*********0**********72
! *** Conservative (flux-form) Semi-Lagrangian transport ***
!****6***0*********0*********0*********0*********0*********0**********72
 
2222  continue

! ghost zone for the western edge:
      iuw =  -c(1,j)
      iuw = min(0, iuw)

      do i=iuw, 0
         qtmp(i) = q(im+i,j)
      enddo

! ghost zone for the eastern edge:
      iue = imp - c(im,j)
      iue = max(imp, iue)

      do i=imp, iue
         qtmp(i) = q(i-im,j)
      enddo

      if(iord.eq.1 .or. j.eq.j1. or. j.eq.j2) then
      do i=1,im
        iu = c(i,j)
      if(c(i,j) .le. 0.) then
        itmp = i - iu
        isave(i) = itmp - 1
      else
        itmp = i - iu - 1
        isave(i) = itmp + 1
      endif
        fx1(i) = (c(i,j)-iu) * qtmp(itmp)
      enddo

! Zero high order contribution
      do i=1,im
         fx2(i,j) = 0.
      enddo

      ELSE
      call xmist(im,IML,qtmp,dc)

      do i=iuw, 0
         dc(i) = dc(im+i)
      enddo

      do i=imp, iue
         dc(i) = dc(i-im)
      enddo

      do i=1,im
            iu  = c(i,j)
            rut = c(i,j) - iu
         if(c(i,j) .le. 0.) then
            itmp = i - iu
            isave(i) = itmp - 1
            fx2(i,j) = -rut*dc(itmp)*(1.+rut)
         else
            itmp = i - iu - 1
            isave(i) = itmp + 1
            fx2(i,j) = rut*dc(itmp)*(1.-rut)
         endif
            fx1(i) = rut*qtmp(itmp)
      enddo

      ENDIF
 
      do i=1,im
      IF(c(i,j).GT.1.) then
!DIR$ NOVECTOR
        do ist = isave(i),i-1
           fx1(i) = fx1(i) + qtmp(ist)
        enddo
      elseIF(c(i,j).LT.-1.) then
!DIR$ NOVECTOR
        do ist = i,isave(i)
           fx1(i) = fx1(i) - qtmp(ist)
        enddo
      endif
      enddo
!DIR$ VECTOR
      do i=1,im
         fx1(i)   = PU(i,j)*fx1(i)
         fx2(i,j) = PU(i,j)*fx2(i,j)
      enddo
 
1309  fx1(IMP  ) = fx1(1  )
      fx2(IMP,j) = fx2(1,j)

! Update using low order fluxes.
      DO i=1,im

         DQ(i,j) =  DQ(i,j) + fx1(i)-fx1(i+1)

! bey, 6/20/00. for mass-flux diagnostic
         fx1_tp(i,j) = fx1(i)
      enddo
 
1310  continue
      return
      end subroutine xtp

!------------------------------------------------------------------------------

      subroutine  ymist(IMR,JNP,j1,P,DC)
!****6***0*********0*********0*********0*********0*********0**********72
      PARAMETER ( R24 = 1./24. )
      REAL P(IMR,JNP),DC(IMR,JNP)
!============================================================================
! Cray NOBOUNDS directive will turn off all subscript bounds checking.
! This directive is also legal on SGI compilers (bmy, 4/24/00)
!DIR$ NOBOUNDS
!============================================================================
!
! 2nd order version for scalars
!
      IMH = IMR / 2
      JMR = JNP - 1
!
      do 10 i=1,IMR*(JMR-1)
      tmp = 0.25*(p(i,3) - p(i,1))
      Pmax = max(p(i,1),p(i,2),p(i,3)) - p(i,2)
      Pmin = p(i,2) - min(p(i,1),p(i,2),p(i,3))
      DC(i,2) = sign(min(abs(tmp),Pmin,Pmax),tmp)
10    CONTINUE
!
! Poles:
!
      if(j1.ne.2) then
      do i=1,IMR
      DC(i,1) = 0.
      DC(i,JNP) = 0.
      enddo
      else
! Determine slopes in polar caps for scalars!
!
      do 20 i=1,IMH
! South
      tmp = 0.25*(p(i,2) - p(i+imh,2))
      Pmax = max(p(i,2),p(i,1), p(i+imh,2)) - p(i,1)
      Pmin = p(i,1) - min(p(i,2),p(i,1), p(i+imh,2))
      DC(i,1)=sign(min(abs(tmp),Pmax,Pmin),tmp)
! North.
      tmp = 0.25*(p(i+imh,JMR) - p(i,JMR))
      Pmax = max(p(i+imh,JMR),p(i,jnp), p(i,JMR)) - p(i,JNP)
      Pmin = p(i,JNP) - min(p(i+imh,JMR),p(i,jnp), p(i,JMR))
      DC(i,JNP) = sign(min(abs(tmp),Pmax,pmin),tmp)
20    continue
!
! Scalars:
      do 25 i=imh+1,IMR
      DC(i,  1) =  - DC(i-imh,  1)
      DC(i,JNP) =  - DC(i-imh,JNP)
25    continue
      endif
      return
      end subroutine ymist

!------------------------------------------------------------------------------

      subroutine ytp(IMR,JNP,j1,j2,acosp,RCAP,DQ,P,C,DC2
     &              ,ymass,fy1,A6,AR,AL,fy2,JORD,
     &              fy1_tp)
!****6***0*********0*********0*********0*********0*********0**********72

      REAL P(IMR,JNP),C(IMR,JNP),ymass(IMR,JNP),fy2(IMR,JNP),
     &     DC2(IMR,JNP),DQ(IMR,JNP),acosp(JNP)

! bey, 6/20/00. for mass-flux diagnostic
      REAL fy1_tp(IMR,JNP)
      
!============================================================================
! Cray NOBOUNDS directive will turn off all subscript bounds checking.
! This directive is also legal on SGI compilers (bmy, 4/24/00)
!DIR$ NOBOUNDS
!============================================================================

! Work array
      REAL fy1(IMR,JNP),AR(IMR,JNP),AL(IMR,JNP),A6(IMR,JNP)
!
      JMR = JNP - 1
      len = IMR*(J2-J1+2)
 
      if(JORD.eq.1) then

      DO 1000 i=1,len
      JT = float(J1) - C(i,J1)
1000  fy1(i,j1) = p(i,JT)

      DO 1050 i=1,len
1050  fy2(i,j1) = 0.

      else
      call ymist(IMR,JNP,j1,P,DC2)
!
      if(JORD.LE.0 .or. JORD.GE.3) then
      call fyppm(C,P,DC2,fy1,fy2,IMR,JNP,j1,j2,A6,AR,AL,JORD)
      else
      DO 1200 i=1,len
      JT = float(J1) - C(i,J1)
      fy1(i,j1) = p(i,JT)
1200  fy2(i,j1) = (sign(1.,C(i,j1))-C(i,j1))*DC2(i,JT)
      endif
      endif
!
      DO 1300 i=1,len
      fy1(i,j1) = fy1(i,j1)*ymass(i,j1)
1300  fy2(i,j1) = fy2(i,j1)*ymass(i,j1)
!
!=============================================================================
! This loop had to be extended for the mass-flux diagnostics (bmy, 4/26/00)
!      DO 1400 j=j1,j2
!      DO 1400 i=1,IMR
!1400  DQ(i,j) = DQ(i,j) + (fy1(i,j) - fy1(i,j+1)) * acosp(j)
!=============================================================================
      DO j=j1,j2
      DO i=1,IMR
         DQ(i,j) = DQ(i,j) + (fy1(i,j) - fy1(i,j+1)) * acosp(j)

         ! bey, 6/20/00. for mass-flux diagnostic
         fy1_tp(i,j) = fy1(i,j) 
#if   defined( TOMAS )
         !Add one more strip of J to the array (win, 7/17/09)
         IF ( j == j2 ) fy1_tp(i,j+1) = fy1(i,j+1)
#endif
      ENDDO
      ENDDO
!
! Poles
      sum1 = fy1(IMR,j1  )
      sum2 = fy1(IMR,J2+1)
      do i=1,IMR-1
      sum1 = sum1 + fy1(i,j1  )
      sum2 = sum2 + fy1(i,J2+1)
      enddo
!
      sum1 = DQ(1,  1) - sum1 * RCAP
      sum2 = DQ(1,JNP) + sum2 * RCAP
      do i=1,IMR
      DQ(i,  1) = sum1
      DQ(i,JNP) = sum2
      enddo
!
      if(j1.ne.2) then
      do i=1,IMR
      DQ(i,  2) = sum1
      DQ(i,JMR) = sum2
      enddo
      endif
      return
      end subroutine ytp

!------------------------------------------------------------------------------

      SUBROUTINE PRESS_FIX( FX, FY, NDT, ACOSP, J1, State_Met )
!
!******************************************************************************
!  Subroutine PRESS_FIX is a wrapper for the Pressure fixer DYN0.  PRESS_FIX 
!  takes the mass fluxes in pressure units and converts them to [kg air/s] 
!  using the correct geometry for TPCORE. (bdf, bmy, 10/11/01, 2/4/03)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) FX     (REAL*8 ) : E-W flux passed from TPCORE [mb/timestep]
!  (2 ) FY     (REAL*8 ) : N-S flux passed from TPCORE [mb/timestep]
!  (3 ) NDT    (INTEGER) : Dynamic timestep for TPCORE [s]
!  (4 ) ACOSP  (REAL*8 ) : Array of inverse cosines    [unitless]
!  (5 ) J1     (INTEGER) : TPCORE polar cap extent     [# of boxes]
! 
!  NOTES:
!  (1 ) Adapted from original code from LLNL.  Added comments and F90 syntax
!        for declarations.  (bdf, bmy, 10/11/01)
!  (2 ) For now, assumes that JJPAR=JJPAR, and DXYP(J) is equivalent to
!        DXYP(J+J0). (bmy, 10/11/01)
!  (3 ) Now declare DXYP as a local array, and initialize it with calls
!        to routine GET_AREA_M2 of "grid_mod.f".  Now use function GET_TS_DYN
!        from "time_mod.f".  Remove reference to CMN header file. (bmy, 2/4/03)
!  01 Mar 2012 - R. Yantosca - Now use GET_AREA_M2(I,J,L) from grid_mod.F90
!******************************************************************************
!
      ! References to F90 modules
      USE GIGC_State_Met_Mod, ONLY : MetState
      USE GRID_MOD,           ONLY : GET_AREA_M2
      USE TIME_MOD,           ONLY : GET_TS_DYN

      USE CMN_SIZE_MOD             ! Size parameters
      USE CMN_DIAG_MOD             ! Diagnostic switches
      USE CMN_GCTM_MOD             ! g0_100
      
      IMPLICIT NONE
#     include "define.h" 

      ! Arguments
      INTEGER,        INTENT(IN)    :: NDT, J1
      REAL*8,         INTENT(IN)    :: ACOSP(JJPAR)
      REAL*8,         INTENT(INOUT) :: FX(IIPAR,JJPAR,LLPAR) 
      REAL*8,         INTENT(INOUT) :: FY(IIPAR,JJPAR,LLPAR)
      TYPE(MetState), INTENT(IN)    :: State_Met   ! Meteorology State object

      ! Local variables
      INTEGER                :: I, J, J2, K, K2, L, AS
      REAL*8                 :: DTC, DTDYN, NSDYN, SUM1, SUM2
      REAL*8                 :: NP_FLUX(IIPAR,LLPAR)
      REAL*8                 :: SP_FLUX(IIPAR,LLPAR)
      REAL*8                 :: ALFA(IIPAR+1,JJPAR,LLPAR)
      REAL*8                 :: BETA(IIPAR,JJPAR+1,LLPAR)
      REAL*8                 :: GAMA(IIPAR,JJPAR,LLPAR+1)
      REAL*8                 :: UMFLX(IIPAR,JJPAR,LLPAR)
      REAL*8                 :: VMFLX(IIPAR,JJPAR,LLPAR)

      ! Local SAVEd variables
!#if defined( DEVEL )
      LOGICAL, SAVE              :: FIRST = .TRUE.
      REAL*8,  ALLOCATABLE, SAVE :: DXYP(:)
!else
!      LOGICAL, SAVE          :: FIRST = .TRUE.
!      REAL*8,  SAVE          :: DXYP(JJPAR)
!#endif
      !=================================================================
      ! PRESS_FIX begins here!
      !
      ! K  is the vertical index down from the atmosphere top downwards
      ! K2 is the vertical index up from the surface
      !=================================================================

      ! Initialize arrays
      ALFA  = 0d0
      BETA  = 0d0
      GAMA  = 0d0

      ! NSDYN is the dynamic time step in seconds
      NSDYN = GET_TS_DYN() * 60d0

      ! J2 is the south polar edge
      J2    = JJPAR - J1 + 1

      ! DTDYN = double precision value for NDT, the dynamic timestep
      DTDYN = DBLE( NDT )

      ! Save grid box surface areas [m2] into the local DXYP array
      IF ( FIRST ) THEN
!#if defined( DEVEL )
         ALLOCATE( DXYP(JJPAR), STAT=AS)
         IF ( AS /=0 ) CALL ALLOC_ERR( 'DXYP_PRESS_FIX')
!#endif

         DO J = 1, JJPAR
            DXYP(J) = GET_AREA_M2( 1, J, 1 )
         ENDDO

         ! Reset first-time flag
         FIRST = .FALSE.
      ENDIF

      !=================================================================
      ! FX is the E-W mass flux from TPCORE in [mb/timestep].  
      ! UMFLX is the mass flux in [kg air/s], which is what DYN0 needs.
      !
      ! FY is the E-W mass flux from TPCORE in [mb/timestep].  
      ! VMFLX is the mass flux in [kg air/s], which is what DYN0 needs.
      !
      ! The unit conversion from [mb/timestep] to [kg air/s] is:
      ! 
      !   mb  | 100 Pa | 1 kg air |  s^2  |  step  | DXYP m^2     kg air
      ! ------+--------+----------+-------+--------+---------- = -------
      !  step |   mb   | Pa m s^2 | 9.8 m | DTDYN s|    s           s
      !=================================================================
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, K, K2, DTC )
      DO K = 1, LLPAR
         K2 = LLPAR - K + 1

         ! Compute UMFLX from FX
         DO J = 1, JJPAR
         DO I = 1, IIPAR
            UMFLX(I,J,K2) = FX(I,J,K) * ( G0_100 * DXYP(J) ) / DTDYN
         ENDDO
         ENDDO

         ! Compute VMFLX from FY
         DO I =  1, IIPAR
         DO J = J1, J2+1
            IF ( FY(I,J,K) .GE. 0 ) THEN
               DTC = FY(I,J,K) * G0_100 * ACOSP(J)  * DXYP(J)   / DTDYN
            ELSE
               DTC = FY(I,J,K) * G0_100 * ACOSP(J-1)* DXYP(J-1) / DTDYN
            ENDIF

            VMFLX(I,J,K2) = DTC
         ENDDO
         ENDDO

         !=================================================================
         ! TREATMENT OF THE POLES: 1
         ! copy ymass values strait into vmflx at poles for pressure fixer
         !=================================================================
         DO I = 1, IIPAR
            VMFLX(I,1,K2)     = FY(I,1,K)
            VMFLX(I,J1-1,K2)  = FY(I,J1-1,K)
            VMFLX(I,JJPAR,K2) = FY(I,JJPAR,K)
         ENDDO
      ENDDO
!$OMP END PARALLEL DO
      
      DO K = 1, LLPAR

         !=================================================================
         ! TREATMENT OF THE POLES: 2
         ! North polar cap: J=1
         !=================================================================
         SUM1 = FY(IIPAR,J1,K)
         DO I = 1, IIPAR-1
            SUM1 = SUM1 + FY(I,J1,K)
         ENDDO

         ! NORTH POLE FLUX IN KG.
         DO I = 1, IIPAR
            NP_FLUX(I,K) = SUM1 * G0_100 * ACOSP(1) * DXYP(1)
         ENDDO 

         !==============================================================
         ! TREATMENT OF THE POLES: 3
         ! South polar cap: J=JJPAR
         !==============================================================
         SUM2 = FY(IIPAR,J2+1,K)
         DO I = 1, IIPAR-1
            SUM2 = SUM2 + FY(I,J2+1,K)
         ENDDO

         DO I = 1, IIPAR
            SP_FLUX(I,K) = SUM2 * G0_100 * ACOSP(JJPAR) * DXYP(JJPAR)
         ENDDO
      ENDDO

      !=================================================================
      ! Call DYN0 to fix the pressures
      !=================================================================
      CALL DYN0( NSDYN, J1,    NP_FLUX, SP_FLUX, 
     &           UMFLX, VMFLX, ALFA,    BETA,   GAMA, State_Met )

      !=================================================================
      ! ALFA is the E-W mass flux adjusted by DYN0 in [kg air/s]
      ! FX   is the E-W mass flux for TPCORE       in [mb/timestep].  
      !
      ! BETA is the N-S mass flux adjusted by DYN0 in [kg air/s]
      ! FY   is the E-W mass flux for TPCORE       in [mb/timestep].  
      !
      ! The unit conversion from to [kg air/s] to [mb/timestep] is:
      ! 
      !  kg air | Pa m s^2 | 9.8 m |    1     | DTDYN s |  mb       mb 
      ! --------+----------+-------+----------+---------+------- = ----
      !    s    | 1 kg air |  s^2  | DXYP m^2 |  step   | 100 Pa   step
      !=================================================================
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, K, K2 )
      DO K = 1, LLPAR
         K2 = LLPAR - K + 1

         ! Update FX from ALFA
         DO J = 1, JJPAR
         DO I = 1, IIPAR
            FX(I,J,K) = ALFA(I,J,K2) * DTDYN / ( G0_100 * DXYP(J) )
         ENDDO
         ENDDO

         ! Update FY from BETA
         DO I =  1, IIPAR
         DO J = J1, J2+1
            IF ( BETA(I,J,K) .GE. 0 ) THEN
               FY(I,J,K) = BETA(I,J,K2) * DTDYN /
     &                     ( G0_100 * ACOSP(J) * DXYP(J) )
            ELSE
                FY(I,J,K) = BETA(I,J,K2) * DTDYN /
     &                     ( G0_100 * ACOSP(J-1) * DXYP(J-1) )
            ENDIF
         ENDDO
         ENDDO

         ! Special treatment of BETA at the poles
         DO I = 1, IIPAR
            FY(I,1,K)     = BETA(I,1,K2)
            FY(I,J1-1,K)  = BETA(I,J1-1,K2)
            FY(I,JJPAR,K) = BETA(I,JJPAR,K2)
         ENDDO
      ENDDO
!$OMP END PARALLEL DO

      ! Return to calling program
      END SUBROUTINE PRESS_FIX

!------------------------------------------------------------------------------

      SUBROUTINE DYN0( DTWIND,   J1,     NP_FLUX,  SP_FLUX,
     &                 UMFLX,    VMFLX,  ALFA,     BETA,     GAMA,
     &                 State_Met )
!
!******************************************************************************
!  Subroutine DYN0 is the pressure fixer for TPCORE.  DYN0 readjusts the
!  mass fluxes ALFA, BETA, GAMA, so that they are consistent with the
!  met fields. (bdf, bmy, 10/11/01, 7/20/04)
!  
!  Arguments as Input:
!  ============================================================================
!  (1 ) DTWIND  (REAL*8 ) : Time step between wind intervals [s]
!  (2 ) J1      (INTEGER) : TPCORE polar cap width 
!  (4 ) NP_FLUX (REAL*8 ) : North polar flux (from PRESS_FIX) in [kg]
!  (5 ) SP_FLUX (REAL*8 ) : South polar flux (from PRESS_FIX) in [kg]
!  (6 ) UMFLX   (REAL*8 ) : Wet air mass flux in E-W     direction [kg air/s]
!  (7 ) VMFLX   (REAL*8 ) : Wet air mass flux in N-S     direction [kg air/s]
!  (8 ) ALFA    (REAL*8 ) : Dry air mass flux in E-W     direction [kg air/s]
!  (9 ) BETA    (REAL*8 ) : Dry air mass flux in N-S     direction [kg air/s]
!  (10) GAMA    (REAL*8 ) : Dry air mass flux in up/down direction [kg air/s]
!
!  Arguments as Output:
!  ============================================================================
!  (8 ) ALFA    (REAL*8 ) : ALFA air mass, after pressure fix is applied
!  (9 ) BETA    (REAL*8 ) : BETA air mass, after pressure fix is applied
!  (10) GAMA    (REAL*8 ) : GAMA air mass, after pressure fix is applied
!
!  NOTES:
!  (1 ) Adapted from original code from LLNL.  Added comments and F90 syntax
!        for declarations.  (bdf, bmy, 10/10/01)
!  (2 ) For a global run (as we usually do in GEOS-CHEM) IM=ID=IIPAR and 
!        JM=JD=JJPAR. (bmy, 10/10/01)
!  (3 ) For now, assumes that JJPAR=JJPAR, and DXYP(J) is equivalent to
!        DXYP(J+J0). (bmy, 10/11/01)
!  (4 ) Rename AD to AD_L so as not to conflict with the AD array in
!        the header file "CMN" (bmy, 10/11/01)
!  (5 ) Now reference PSC2 instead of PS from "dao_mod.f".  Replace all 
!        instances of PS with PSC2.  Updated comments. (bdf, bmy, 4/1/02)
!  (6 ) Removed obsolete code from 4/1/02 (bdf, bmy, 8/22/02)
!  (7 ) Now declare DXYP as a local array, and initialize it with calls
!        to routine GET_AREA_M2 and GET_YOFFSET of "grid_mod.f".  Now also
!        references GET_BP from "pressure_mod.f" (bmy, 2/11/03)
!  (8 ) Removed reference to CMN (bmy, 7/20/04)
!  01 Mar 2012 - R. Yantosca - Now use GET_AREA_CM2(I,J,L) from grid_mod.F90
!  09 Nov 2012 - M. Payer    - Replaced all met field arrays with State_Met
!                              derived type object
!******************************************************************************
!
      ! References to F90 modules
      USE GIGC_State_Met_Mod, ONLY : MetState
      USE GRID_MOD,           ONLY : GET_AREA_M2, GET_YOFFSET
      USE PRESSURE_MOD,       ONLY : GET_BP

      USE CMN_SIZE_MOD             ! Size parameters

      IMPLICIT NONE
#     include "define.h"

      ! Arguments
      INTEGER,        INTENT(IN)    :: J1 
      REAL*8,         INTENT(IN)    :: DTWIND
      REAL*8,         INTENT(IN)    :: NP_FLUX(IIPAR,LLPAR)
      REAL*8,         INTENT(IN)    :: SP_FLUX(IIPAR,LLPAR)
      REAL*8,         INTENT(IN)    :: UMFLX(IIPAR,JJPAR,LLPAR)
      REAL*8,         INTENT(IN)    :: VMFLX(IIPAR,JJPAR,LLPAR)
      REAL*8,         INTENT(INOUT) :: ALFA(IIPAR+1,JJPAR,LLPAR)
      REAL*8,         INTENT(INOUT) :: BETA(IIPAR,JJPAR+1,LLPAR)
      REAL*8,         INTENT(INOUT) :: GAMA(IIPAR,JJPAR,LLPAR+1)
      TYPE(MetState), INTENT(IN)    :: State_Met   ! Meteorology State object

      ! Local variables
      LOGICAL                :: LSP, LNP, LEW
      INTEGER                :: IIX, JJX, KM, JB, JE, IEPZ, IMZ
      INTEGER                :: I, J, J2, K, L, AS
      REAL*8                 :: ALFAX, UFILT, VFILT, PCTM8, G100
      REAL*8                 :: AIRQAV, AWE, SUMAD0, SUMAW0, AIRWET
      REAL*8                 :: AIRH2O, AIRQKG ,SUM1, SUMA, SUMP, SUMQ
      REAL*8                 :: SUMU, SUMV, SUMW, ZIMZ, ZDTW, G0
      REAL*8                 :: AD_L(IIPAR,JJPAR,LLPAR) 
      REAL*8                 :: AIRD(IIPAR,JJPAR,LLPAR)
      REAL*8                 :: AIRNEW(IIPAR,JJPAR,LLPAR)
      REAL*8                 :: AIRX(IIPAR,JJPAR,LLPAR)
      REAL*8                 :: AX(IIPAR+1,JJPAR)
      REAL*8                 :: BX(IIPAR,JJPAR+1)
      REAL*8                 :: MERR(IIPAR,JJPAR)
      REAL*8                 :: PCTM(IIPAR,JJPAR)
      REAL*8                 :: PERR(IIPAR,JJPAR)
      REAL*8                 :: SPHU_KG(IIPAR,JJPAR,LLPAR)
      REAL*8                 :: SUMAQ(IIPAR,JJPAR)
      REAL*8                 :: XYB(IIPAR,JJPAR)
      REAL*8                 :: XYZB(IIPAR,JJPAR,LLPAR)

      ! Local saved variables
      LOGICAL, SAVE          :: FIRST = .TRUE.

!#if defined( DEVEL )
      REAL*8,  ALLOCATABLE, SAVE :: DXYP(:)
      REAL*8,  ALLOCATABLE, SAVE :: DSIG(:)
!#else
!      REAL*8,  SAVE          :: DXYP(JJPAR)
!      REAL*8,  SAVE          :: DSIG(LLPAR)
!#endif

      !=================================================================
      ! DYN0 begins here!
      !
      ! UNITS OF AIR MASS AND TRACER = (kg)  
      !
      ! Air mass (kg) is given by:
      !    area (m^2) * pressure thickness (Pa) / g0
      !
      ! DXYP(J)        = area of [I,J] [m^2] (is longitude-symmetric)
      !
      ! PSC2(I,J)      = surf pressure [Pa] averaged in extended zone.
      !                  this is the surface pressure at the end of the
      !                  current dynamic timestep, passed to TPCORE.
      !
      ! SPHU_KG(I,J,K) = specific humidity of grid box 
      !                  [kg H2O/kg wet air] averaged in extended zone.
      !
      ! AIRQKG(I,J)    = Mass of H2O [kg] at each level
      !                = PS(I,J)) * SPHU_KG(I,J,K)
      !
      ! AIRD(I,J,K)    = dry-air mass [kg] in each box as calculated 
      !                  in CTM at the beginning of each time step, 
      !                  updated at end of DYN0.
      !
      ! PCTM(I,J)      = inferred wet-air (total) surf press [Pa] calc. 
      !                  in CTM (using SUMAQ & AIRD-X-NEW)
      !
      ! AIRNEW(I,J,K)  = new dry-air mass in each CTM box after 
      !                  horizontal divergence (ALFA+BETA) over time 
      !                  step DTWIND (sec)
      !
      ! AIRX(I,J,K)    = expected dry-air mass in each CTM box after 
      !                  calculating the vertical divergence (GAMA)   
      !                  (also used for GCM dry mass)
      !                = XYZA(I,J,K) + XYZB(I,J,K)*PCTM(I,J) - AIRQKG
      !
      ! DTWIND         = time step [s] that applies to the averaged 
      !                  wind fields (i.e., the time between successive 
      !                  pressures.
      !
      !-----------------------------------------------------------------
      !
      ! Assume that we have "wet-air" mass fluxes across each boundary
      !
      !   UMFLX(I,J,K)   ==> [I,J,K] ==>  UMFLX(I+1,J,K)   [kg air/s]
      !   VMFLX(I,J,K)   ==> [I,J,K] ==>  VMFLX(I,J+1,K)   [kg air/s]
      ! 
      ! Convert to "dry-air" mass flux in/out of box using 
      ! average Q at boundary
      !
      !   ALFA(I,J,K)   ==> [I,J,K] ==>  ALFA(I+1,J,K)     [kg air/s]
      !   BETA(I,J,K)   ==> [I,J,K] ==>  BETA(I,J+1,K)     [kg air/s]
      ! 
      ! Calculate convergence in each layer of dry air, compare with 
      ! expected dry air mass (AIRX) and then calculate vertical 
      ! dry-mass fluxes
      !
      !   GAMA(I,J,K)   ==> [I,J,K] ==>  GAMA(I,J,K+1)     [kg air/s]
      ! 
      ! Horizontal pressure filter adjusts UMFLX & VMFLX to reduce 
      ! error in [PCTM - PSC2] 
      !
      !     UMFLX + pressure filter ==> UMFLX#,  
      !     VMFLX + filter ==> VMFLX# (temporary)
      !
      ! The pressure filter does nearest neighbor flux 
      ! (adjusting ALFA/BETA)
      !
      !-----------------------------------------------------------------
      ! 
      ! Note that K->K+1 is downward (increasing pressure) and 
      ! that boundaries:
      !   GAMA(I,J,1) = GAMA(I,J,KM+1)    = 0   no flux across 
      !                                         upper/lower boundaries
      !
      !   BETA(I,1,K) = BETA(I,JJPAR+1,K) = 0   no flux at S & N poles
      !
      !   ALFA(1,J,K) = ALFA(IIPAR+1,J,K)       is NOT ZERO, but cyclic
      !
      ! Dimensions for ALFA, BETA, GAMA are extended by +1 beyond grid 
      ! to allow simple formulation of fluxes in/out of final grid box.
      ! 
      ! GCM input UMFLX,VMFLX,PSG is ALWAYS of GLOBAL dimensions 
      ! (IIPAR x JJPAR x LLPAR)
      !
      ! Indices of ALFA, BETA, GAMA, SPHU_KG & PSC2 are always LOCAL 
      ! (IIPAR x JJPAR x KM): FOR GEOS-CHEM, KM = LLPAR (bmy
      !
      ! Indices of tracer (STT), and diagnostics are local 
      ! (w.r.t. WINDOW.  WINDOW calculations are defined by an 
      ! offset and size
      !
      !         I0 .ge.0 and IIPAR+I0 .le. IIPAR
      !         J0 .ge.0 and JJPAR+J0 .le. JJPAR
      !         K0 .ge.0 and KM+K0 .le. LLPAR
      !
      ! The WINDOW calculation must allow for a boundary layer 
      ! of grid boxes:
      !
      !         IG(abs. coords) = IW(in window) + I0
      !         JG(abs. coords) = JW(in window) + J0
      !         KG(abs. coords) = KW(in window) + K0
      !
      ! vertical window (NEW) allows for an upper boundary with flow 
      ! across it and specified mixing ratio b.c.'s at KG = K0
      !=================================================================

      ! First-time initialization
      IF ( FIRST ) THEN
         
         ! Surface area [m2]
#if defined( DEVEL )
         ALLOCATE( DXYP(JJPAR), STAT=AS)
         IF ( AS /=0 ) CALL ALLOC_ERR( 'DXYP_DYN0')
#endif

         DO J = 1, JJPAR
            DXYP(J) = GET_AREA_M2( 1, J, 1 )
         ENDDO

         ! Sigma-level thickness [unitless]
         ! Assumes we are using a pure-sigma grid
#if defined( DEVEL )
         ALLOCATE( DSIG(LLPAR), STAT=AS)
         IF ( AS /=0 ) CALL ALLOC_ERR( 'DSIG')
#endif

         DO L = 1, LLPAR        
            DSIG(L) = GET_BP(L) - GET_BP(L+1)
         ENDDO

         ! Reset first-time flag
         FIRST = .FALSE.
      ENDIF

      ! for tpcore poles.
      J2 = JJPAR - J1 + 1

      ! geos code
      G0 = 9.8d0 

      !=================================================================
      ! XYZB is the factor needed to get mass in kg of gridbox
      ! mass (kg) = XYZB (kg/mb) * P (mb)
      !
      ! AD_L is the dry air mass in the grid box
      !
      ! SPHU_KG is the water vapor [kg H2O/kg air]
      !=================================================================
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L )
      DO L = 1, LLPAR
      DO J = 1, JJPAR
      DO I = 1, IIPAR
         XYZB(I,J,L)    = DSIG(L) * DXYP(J) * 1.d2 / G0
         AD_L(I,J,L)    = State_Met%AIRDEN(L,I,J) *
     &                    State_Met%AIRVOL(I,J,L)
         SPHU_KG(I,J,L) = State_Met%SPHU(I,J,L) / 1000d0
      ENDDO
      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      !=================================================================
      ! XYB is the factor needed to get mass in kg of column
      !=================================================================
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J )
      DO J = 1, JJPAR
      DO I = 1, IIPAR
         XYB(I,J) = SUM( XYZB(I,J,1:LLPAR) )
      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      !=================================================================
      ! Define other variables
      !=================================================================
      G100   = 100.D0 / G0
      ZDTW   =   1.D0 / DTWIND
      LSP    = (         GET_YOFFSET() .EQ. 0     )
      LNP    = ( JJPAR + GET_YOFFSET() .EQ. JJPAR )
      LEW    = (                 IIPAR .EQ. IIPAR )

      !=================================================================
      ! Initialize ALFA with UMFLX and BETA with VMFLX
      !=================================================================
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L )
      DO L = 1, LLPAR
         DO J = 1, JJPAR
            DO I = 1,IIPAR
               ALFA(I,J,L) = UMFLX(I,J,L)
            ENDDO

            ALFA(IIPAR+1,J,L) = ALFA(1,J,L)
         ENDDO

         DO J = 2, JJPAR
            DO I = 1, IIPAR
               BETA(I,J,L) = VMFLX(I,J,L)
            ENDDO

            DO I = 1, IIPAR
               BETA(I,1,L)       = 0.D0
               BETA(I,JJPAR+1,L) = 0.D0
            ENDDO
         ENDDO
      ENDDO
!$OMP END PARALLEL DO

      !=================================================================
      ! SUMAQ(I,J): column integral of water (kg)
      ! Check on air mass
      !=================================================================
      SUMAD0 = 0.D0
      SUMAW0 = 0.D0

      DO J = 1, JJPAR
      DO I = 1, IIPAR
         SUMAQ(I,J) = 0.D0
        
         DO K = 1, LLPAR
            AIRWET     = State_Met%PSC2(I,J) * XYZB(I,J,K)
            AIRH2O     = SPHU_KG(I,J,K) * AIRWET
            SUMAQ(I,J) = SUMAQ(I,J)     + AIRH2O
            SUMAD0     = SUMAD0         + AIRWET
            SUMAW0     = SUMAW0         + AIRH2O
         ENDDO
      ENDDO
      ENDDO

      SUMAD0 = SUMAD0 - SUMAW0

      !=================================================================
      ! Initialize AIRD, the dry-air mass [kg] in each box as calculated 
      ! in CTM at the start of each time step, updated at end of DYN0.
      !
      ! Compute AIRNEW, the new dry-air mass in each CTM box after 
      ! horizontal divergence (ALFA+BETA) over time step DTWIND (sec)
      !=================================================================
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, K )
      DO K =  1, LLPAR
      DO J = J1, J2
      DO I =  1, IIPAR
         AIRD(I,J,K)   = AD_L(I,J,K)
         AIRNEW(I,J,K) = AIRD(I,J,K) + DTWIND *
     &                                 ( ALFA(I,J,K) - ALFA(I+1,J,K) +
     &                                   BETA(I,J,K) - BETA(I,J+1,K) )
      ENDDO
      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      !=================================================================
      ! treatment of the poles for tpcore.
      ! j=2 and j=jjpar-1 don't have any airmass change.
      !=================================================================
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, K )
      DO K = 1, LLPAR

         ! J=1
         DO I = 1, IIPAR
            AIRNEW(I,1,K) = AD_L(I,1,K) - NP_FLUX(I,K)
         ENDDO
       
         ! J=JJPAR
         DO I = 1, IIPAR
            AIRNEW(I,JJPAR,K) = AD_L(I,JJPAR,K) + SP_FLUX(I,K)
         ENDDO
      ENDDO
!$OMP END PARALLEL DO

      !=================================================================
      ! Average AIRNEW at the South pole
      !=================================================================
      ZIMZ = 1.D0 / DBLE( IIPAR )

      IF ( LSP ) THEN
        JB = 2

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, K, SUMA )
        DO K = 1, LLPAR
           SUMA = SUM( AIRNEW(1:IIPAR,1,K) ) * ZIMZ

           DO I = 1, IIPAR
              AIRNEW(I,1,K) = SUMA
           ENDDO
        ENDDO
!$OMP END PARALLEL DO
      ELSE
         JB = 1
      ENDIF

      !=================================================================
      ! Average AIRNEW at the North pole
      !=================================================================
      IF ( LNP ) THEN
        JE = JJPAR - 1

        ! poles, just average AIRNEW
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, K, SUMA )
        DO K = 1, LLPAR
           SUMA = SUM( AIRNEW(1:IIPAR,JJPAR,K) ) * ZIMZ

           DO I = 1, IIPAR
              AIRNEW(I,JJPAR,K) = SUMA
           ENDDO
        ENDDO
!$OMP END PARALLEL DO
      ELSE
         JE = JJPAR
      ENDIF

      !================================================================
      ! BEGIN FILTER of PRESSURE ERRORS
      !
      ! Define the error in surface pressure PERR expected at end of 
      ! time step filter by error in adjacent boxes, weight by areas, 
      ! adjust ALFA & BETA
      !
      ! PCTM(I,J) = new CTM wet-air column based on 
      !             dry-air convergence (Pascals) 
      ! PERR(I,J) = pressure-error between CTM-GCM at new time 
      !             (before filter)
      !================================================================
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J )
      DO J = 1, JJPAR
      DO I = 1, IIPAR
         PCTM(I,J) = SUM( AIRNEW(I,J,:) ) / XYB(I,J)
         
         ! special case for j=2, jjpar-1 for tpcore pole configuration.
         IF ( J .eq. 2 .OR. J .eq. JJPAR-1 ) THEN
            PCTM(I,J) = State_Met%PSC2(I,J)
         ENDIF

         PERR(I,J) = PCTM(I,J) - State_Met%PSC2(I,J)
         MERR(I,J) = PERR(I,J) * DXYP(J) * G100
      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      ! Call pressure filter
      CALL PFILTR( MERR,  AX,    BX, DXYP, IIPAR, JJPAR,
     &             IIPAR, JJPAR, 1,  LSP,  LNP,   LEW )

      !=================================================================
      ! Calculate corrections to ALFA from the filtered AX
      !=================================================================
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, IIX, J, K, UFILT )
      DO J = JB, JE
      DO I =  1, IIPAR+1
         IIX   = MIN(I,IIPAR)
         UFILT = AX(I,J) / ( XYB(IIX,J) * DTWIND )

         DO K = 1, LLPAR
            ALFA(I,J,K) = ALFA(I,J,K) + UFILT * XYZB(IIX,J,K)
         ENDDO
      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      !=================================================================
      ! Calculate corrections to BETA from the filtered BX
      !=================================================================
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, JJX, K, VFILT )
      DO J = 1, JJPAR+1
         JJX = J
         IF ( J+J .gt. JJPAR )  JJX = J - 1

         DO I = 1, IIPAR
            VFILT = BX(I,J) / ( XYB(I,JJX) * DTWIND )

            DO K = 1, LLPAR
               BETA(I,J,K) = BETA(I,J,K) + VFILT * XYZB(I,JJX,K)
            ENDDO
         ENDDO
      ENDDO
!$OMP END PARALLEL DO

      !=================================================================
      ! Calculate the corrected AIRNEW's & PCTM after P-filter:
      ! has changed ALFA+BETAs and ctm surface pressure (PCTM)
      !=================================================================
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, K )
      DO K = 1, LLPAR
      DO J = 1, JJPAR
      DO I = 1, IIPAR
         AIRNEW(I,J,K) = AIRD(I,J,K) + DTWIND *
     &                                 ( ALFA(I,J,K) - ALFA(I+1,J,K) +
     &                                   BETA(I,J,K) - BETA(I,J+1,K) )
      ENDDO
      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      !=================================================================
      ! Average the adjusted AIRNEW at the South pole
      !=================================================================      
      ZIMZ = 1.D0 / DBLE( IIPAR ) 

      IF ( LSP ) THEN
         JB = 2
         
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, K, SUMA )
         DO K = 1, LLPAR
            SUMA = SUM( AIRNEW(1:IIPAR,1,K ) ) * ZIMZ 

            DO I = 1, IIPAR
               AIRNEW(I,1,K) = SUMA
            ENDDO
         ENDDO
!$OMP END PARALLEL DO
      ELSE
         JB = 1
      ENDIF

      !=================================================================
      ! Average the adjusted AIRNEW at the North pole
      !=================================================================     
      IF ( LNP ) THEN
         JE = JJPAR -1

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, K, SUMA )
         DO K = 1, LLPAR
            SUMA = SUM( AIRNEW(1:IIPAR,JJPAR,K) ) * ZIMZ

            DO I = 1,IIPAR
               AIRNEW(I,JJPAR,K) = SUMA
            ENDDO
         ENDDO
!$OMP END PARALLEL DO
      ELSE
         JE = JJPAR
      ENDIF
      
      !=================================================================
      ! END OF PRESSURE FILTER
      !
      ! GAMA:  redistribute the new dry-air mass consistent with the
      ! new CTM surface pressure, rigid upper b.c., no change in PCTM
      !
      ! AIRX(I,J,K) = dry-air mass expected, based on PCTM
      !               PCTM(I,J) & PERR(I,J)
      !=================================================================

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, K, PCTM8, AIRQKG )
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         PCTM8     = ( SUM( AIRNEW(I,J,:) ) + SUMAQ(I,J) ) / XYB(I,J)
         PCTM(I,J) = PCTM8
         PERR(I,J) = PCTM8 - State_Met%PSC2(I,J)

         DO K = 1, LLPAR
            AIRQKG      = SPHU_KG(I,J,K) *
     &                  ( XYZB(I,J,K) * State_Met%PSC2(I,J) )
            AIRX(I,J,K) = PCTM8 * XYZB(I,J,K) - AIRQKG
         ENDDO
      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      !=================================================================
      ! GAMA from top down to be consistent with AIRX, AIRNEW not reset!
      !=================================================================
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, K )
      DO J = 1, JJPAR
      DO I = 1, IIPAR
         GAMA(I,J,LLPAR+1) = 0.D0

         DO K = LLPAR, 2, -1
            GAMA(I,J,K) = GAMA(I,J,K+1) - (AIRNEW(I,J,K) - AIRX(I,J,K))
         ENDDO

         ! GAMA(I,J,1) will not be exactly ZERO, but it must be set so!
         GAMA(I,J,1) = 0.D0

         DO K = 2, LLPAR
            GAMA(I,J,K) = GAMA(I,J,K) * ZDTW
         ENDDO
      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      ! Return to calling program
      END SUBROUTINE DYN0

!------------------------------------------------------------------------------

      SUBROUTINE PFILTR( MERR, ALFAX, BETAX, AXY, ID,  JD, 
     &                   IM,   JM,    NITR,  LSP, LNP, LEW )
!
!******************************************************************************
!  Subroutine PFILTR applies the pressure-filter, the pressure 
!  between predicted Ps(CTM) and Ps(GCM). (bdf, bmy, 10/11/01)
!  
!  Arguments as Input:
!  ============================================================================
!  (1  ) MERR(ID,JD)    (REAL*8 ) : mass error
!  (2  ) ALFAX(ID+1,JD) (REAL*8 ) : perturbed ALFA by MERR
!  (3  ) BETAX(ID,JD+1) (REAL*8 ) : perturbed BETA by MERR
!  (4  ) AXY(ID,JD)     (REAL*8 ) : area of grid box (I,J) in [m^2]
!  (5-6) ID, JD         (INTEGER) : "Global" array dimensions for lon, lat
!  (7-8) IM, JM         (INTEGER) : "Window" array dimensions for lon, lat
!  (9  ) NITR           (INTEGER) : number of iterations (NITR .LE. 4)
!  (10 ) LSP            (LOGICAL) : true if J=1  is S. POLE
!  (11 ) LNP            (LOGICAL) : true if J=JM is N. POLE
!  (12 ) LEW            (LOGICAL) : true if cyclic in W-E direction
!                                   (i.e. if I=1 connects to I=IM)
!
!  Arguments as Output:
!  ============================================================================
!  (1  ) MERR(ID,JD)    (REAL*8 ) : adjusted mass error
!  (2  ) ALFAX(ID+1,JD) (REAL*8 ) : adjusted ALFAX 
!  (3  ) BETAX(ID,JD+1) (REAL*8 ) : adjusted BETAX 
!
!  NOTES:
!  (1 ) Adapted from original code from LLNL.  Added comments and F90 syntax
!        for declarations. (bdf, bmy, 10/1/01)
!  (2 ) For a global run (as we usually do in GEOS-CHEM) IM=ID=IIPAR and 
!        JM=JD=JJPAR. (bmy, 10/11/01)
!  (3 ) Removed IMEPZ -- we don't need this for GEOS-CHEM. (bmy, 10/18/01)
!******************************************************************************
!
      IMPLICIT NONE
#     include "define.h"

      ! Arguments
      LOGICAL, INTENT(IN)    :: LSP,LNP,LEW
      INTEGER, INTENT(IN)    :: ID, JD, IM, JM, NITR
      REAL*8,  INTENT(IN)    :: AXY(JD)
      REAL*8,  INTENT(INOUT) :: MERR(ID,JD)
      REAL*8,  INTENT(INOUT) :: ALFAX(ID+1,JD)
      REAL*8,  INTENT(INOUT) :: BETAX(ID,JD+1)

      ! Local variables
      LOGICAL                :: LPOLE
      INTEGER                :: I, J, K 
      REAL*8                 :: X0(ID,JD)

      !=================================================================
      ! PFILTR begins here!
      !=================================================================

      ! LPOLE is true if J=1 is the SOUTH POLE and J=JM is the NORTH POLE
      ! (this is the way GEOS-CHEM is set up, so LPOLE should be TRUE!)
      LPOLE = ( LSP .AND. LNP )

      ! Zero ALFAX, BETAX, save MERR in X0
      DO J = 1, JM
         DO I = 1, IM
            ALFAX(I,J) = 0.D0
            BETAX(I,J) = 0.D0
            X0(I,J)    = MERR(I,J) 
         ENDDO

         ALFAX(IM+1,J) = 0.D0
      ENDDO

      DO I = 1, IM
         BETAX(I,JM+1) = 0.D0
      ENDDO

      !=================================================================
      ! Call LOCFLT to do the local filtering
      !=================================================================
      CALL LOCFLT( MERR, ALFAX, BETAX, AXY, ID,  JD,
     &             IM,   JM,    5,     LSP, LNP, LEW )

      !=================================================================
      ! Call POLFLT to do the pole filtering (if necessary)
      !=================================================================
      IF ( LPOLE ) THEN
         CALL POLFLT( MERR, BETAX, AXY, 1.D0, ID, JD, IM, JM )
      ENDIF
      
      !=================================================================
      ! Compute mass error MERR and return
      ! MERR, ALFAX, and BETAX are now adjusted 
      !=================================================================
      DO J = 1, JM
      DO I = 1, IM
         MERR(I,J) = X0(I,J) + ALFAX(I,J) - ALFAX(I+1,J)
     &                       + BETAX(I,J) - BETAX(I,J+1)
      ENDDO
      ENDDO

      ! Return to calling program
      END SUBROUTINE PFILTR

!------------------------------------------------------------------------------
      
      SUBROUTINE LOCFLT( XERR, AX, BX,   AXY, ID,  JD, 
     &                   IM,   JM, NITR, LSP, LNP, LEW )
!
!******************************************************************************
!  Subroutine LOCFLT applies the pressure-filter to non-polar boxes.  
!  LOCFLT is called from subroutine PFILTR above (bdf, bmy, 10/11/01)
!  
!  Arguments as Input:
!  ============================================================================
!  (1  ) XERR(ID,JD) (REAL*8 ) : mass error
!  (2  ) AX(ID+1,JD) (REAL*8 ) : perturbed ALFA by XERR
!  (3  ) BX(ID,JD+1) (REAL*8 ) : perturbed BETA by XERR
!  (4  ) AXY(ID,JD)  (REAL*8 ) : area of grid box (I,J) in [m^2]
!  (5-6) ID, JD      (INTEGER) : "Global" array dimensions for lon, lat
!  (7-8) IM, JM      (INTEGER) : "Window" array dimensions for lon, lat
!  (9  ) NITR        (INTEGER) : number of iterations (NITR .LE. 4)
!  (10 ) LSP         (LOGICAL) : true if J=1  is S. POLE
!  (11 ) LNP         (LOGICAL) : true if J=JM is N. POLE
!  (12 ) LEW         (LOGICAL) : true if cyclic in W-E direction 
!                                (i.e. if I=1 connects to I=IM)
!
!  Arguments as Output:
!  ============================================================================
!  (1  ) XERR(ID,JD) (REAL*8 ) : adjusted mass error
!  (2  ) AX(ID+1,JD) (REAL*8 ) : adjusted AX 
!  (3  ) BX(ID,JD+1) (REAL*8 ) : adjusted BX
!
!  NOTES:
!  (1 ) Adapted from original code from LLNL.  Added comments and F90 syntax
!        for declarations. (bdf, bmy, 10/11/01)
!  (2 ) For a global run (as we usually do in GEOS-CHEM) IM=ID=IIPAR and 
!        JM=JD=JJPAR. (bmy, 10/11/01)
!******************************************************************************
!
      IMPLICIT NONE
#     include "define.h"
      
      ! Arguments
      LOGICAL, INTENT(IN)    :: LSP, LNP, LEW
      INTEGER, INTENT(IN)    :: ID, JD, IM, JM, NITR
      REAL*8,  INTENT(IN)    :: AXY(JD)
      REAL*8,  INTENT(INOUT) :: XERR(ID,JD)
      REAL*8,  INTENT(INOUT) :: AX(ID+1,JD)
      REAL*8,  INTENT(INOUT) :: BX(ID,JD+1)

      ! Local variables
      INTEGER                :: I, IA, NAZ, J, J1, J2, NFLTR
      REAL*8                 :: SUMA, FNAZ8
      REAL*8                 :: X0(ID,JD)

      !=================================================================
      ! LOCFLT begins here!
      !
      ! Initialize corrective column mass flows (kg):  AX->alfa, BX->beta
      !=================================================================
      DO J = 1, JM
      DO I = 1, IM
         X0(I,J) = XERR(I,J)
      ENDDO
      ENDDO

      !=================================================================
      ! Iterate over mass-error filter
      ! accumulate corrections in AX & BX
      !=================================================================
      DO NFLTR = 1, NITR

         !==============================================================
         ! calculate AX = E-W filter
         !==============================================================

         ! Compute polar box limits
         J1 = 1
         J2 = JM
         IF ( LSP ) J1 = 2
         IF ( LNP ) J2 = JM - 1

         ! Loop over non-polar latitudes
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, FNAZ8 )
         DO J = J1, J2

            ! Calculate pressure-filter E-W wind between boxes [I-1] & [I].  
            ! Enhance filtered wind by size of EPZ, will redistribute 
            ! later within
            FNAZ8 = 0.125d0

            DO I = 2, IM
               AX(I,J) = AX(I,J) + FNAZ8 *(XERR(I-1,J) - XERR(I,J))
            ENDDO

            ! calculate pressure-filter E-W wind at edges I=1 & I=IM+1
            IF ( LEW )  THEN
               AX(IM+1,J) = AX(IM+1,J) + FNAZ8 * (XERR(IM,J) -XERR(1,J))
               AX(1,J)    = AX(1,J)    + FNAZ8 * (XERR(IM,J) -XERR(1,J))
            ELSE                
               ! WINDOW, assume zero error outside window
               AX(1,J)   = AX(1,J)    - FNAZ8 * XERR(1,J)
               AX(IM+1,J)= AX(IM+1,J) + FNAZ8 * XERR(IM,J)
            ENDIF
         ENDDO
!$OMP END PARALLEL DO

         !==============================================================
         ! calculate BX = N-S filter, N-S wind between boxes [J-1] & [J]
         !==============================================================
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, FNAZ8 )
         DO J = 3, JM-1
            FNAZ8 = 0.25D0 * AXY(J) / ( AXY(J-1) + AXY(J) )

            DO I = 1, IM
               BX(I,J) = BX(I,J) + FNAZ8 * ( XERR(I,J-1) - XERR(I,J) )
            ENDDO
         ENDDO
!$OMP END PARALLEL DO

         ! enhance the filtering by factor of 2 ONLY into/out-of polar caps
         FNAZ8 = 0.5D0 * AXY(2) / ( AXY(1) + AXY(2) )

         ! When LSP=TRUE then J=1 is SOUTH POLE
         IF ( LSP ) THEN
            DO I = 1, IM
               BX(I,2) = BX(I,2) + FNAZ8 * (XERR(I,1) -XERR(I,2))
            ENDDO
         ELSE
            DO I = 1, IM
               BX(I,1)= BX(I,1) -0.5D0 *FNAZ8 * XERR(I,1)
               BX(I,2)= BX(I,2) +0.5D0 *FNAZ8 * (XERR(I,1) - XERR(I,2))
            ENDDO
         ENDIF

         FNAZ8 = 0.5D0 * AXY(JM) / ( AXY(JM-1) + AXY(JM) )

         ! When LNP=TRUE, then J=JM is NORTH POLE
         IF ( LNP )  THEN
            DO I = 1, IM
               BX(I,JM) = BX(I,JM) +FNAZ8 *(XERR(I,JM-1) -XERR(I,JM))
            ENDDO
         ELSE
            DO I = 1,IM
               BX(I,JM+1)= BX(I,JM+1) + 0.5D0 *FNAZ8 * XERR(I,JM)
               BX(I,JM)  = BX(I,JM)   + 0.5D0 *FNAZ8 * 
     &                                  (XERR(I,JM-1) -XERR(I,JM))
            ENDDO
         ENDIF

         !==============================================================
         ! need N-S flux across boundaries if window calculation
         ! (assume XERR=0 outside)
         !
         ! JM for optimal matrix/looping, it would be best to 
         ! define XERR=0 for an oversized array XERR(0:IM+1,0:JM+1)
         ! Update the mass error (XERR)
         !==============================================================
         DO J = 1, JM
         DO I = 1, IM
            XERR(I,J) = X0(I,J) + AX(I,J) - AX(I+1,J) 
     &                          + BX(I,J) - BX(I,J+1)
         ENDDO
         ENDDO

      ENDDO  ! NFLTR

      ! Return to calling program
      END SUBROUTINE LOCFLT

!------------------------------------------------------------------------------

      SUBROUTINE POLFLT( XERR, BX, AXY, COEF, ID, JD, IM, JM )
!
!******************************************************************************
!  Subroutine POLFLT applies the pressure-filter to polar boxes.  
!  POLFLT is called from subroutine PFILTR above (bdf, bmy, 10/10/01)
!  
!  Arguments as Input:
!  ============================================================================
!  (1  ) XERR(ID,JD) (REAL*8 ) : mass error
!  (2  ) BX(ID,JD+1) (REAL*8 ) : perturbed BETA by XERR
!  (3  ) AXY(ID,JD)  (REAL*8 ) : area of grid box (I,J) in [m^2]
!  (4  ) COEF        (REAL*8 ) : Multiplicative coefficient ?????
!  (5-6) ID, JD      (INTEGER) : "Window" array dimensions for lon, lat
!  (7-8) IM, JM      (INTEGER) : "Global" array dimensions for lon, lat
!
!  Arguments as Output:
!  ============================================================================
!  (1  ) XERR(ID,JD) (REAL*8 ) : adjusted mass error
!  (2  ) BX(ID,JD+1) (REAL*8 ) : adjusted BX
!
!  NOTES:
!  (1 ) Adapted from original code from LLNL.  Added comments and F90 syntax
!        for declarations. (bdf, bmy, 10/10/01)
!  (2 ) For a global run (as we usually do in GEOS-CHEM) IM=ID=IIPAR and 
!        JM=JD=JJPAR. (bmy, 10/20/01)
!******************************************************************************
!
      IMPLICIT NONE
#     include "define.h"

      ! Arguments
      INTEGER, INTENT(IN)    :: ID, JD, IM, JM
      REAL*8,  INTENT(IN)    :: AXY(JD)
      REAL*8,  INTENT(IN)    :: COEF
      REAL*8,  INTENT(INOUT) :: XERR(ID,JD)
      REAL*8,  INTENT(INOUT) :: BX(ID,JD+1)

      ! Local variables
      INTEGER                :: I, J
      REAL*8                 :: ERAV, BXJ(JD+1), TOTAL
      
      !================================================================= 
      ! POLFLT begins here!
      !
      ! Initialize corrective column mass flows (kg):  BXJ->beta
      !================================================================= 
      DO I = 1, IM

         ! Initialize
         ERAV  = 0.D0
         TOTAL = 0.D0

         ! Sum XERR in ERAV and sum AXY in TOTAL
         DO J = 1, JM
            ERAV  = ERAV  + XERR(I,J)
            TOTAL = TOTAL + AXY(J)
         ENDDO

         ! Compute area-weighted mass error total
         ERAV = ERAV / TOTAL

         ! mass-error filter, make corrections in BX
         BXJ(1) = 0.D0

         DO J = 2, JM
            BXJ(J) = BXJ(J-1) + XERR(I,J-1) - AXY(J-1) * ERAV
         ENDDO

         DO J = 2, JM
            BX(I,J) = BX(I,J) + COEF * BXJ(J)
         ENDDO
      
      ENDDO  ! I

      ! Update XERR
      DO J = 1, JM
      DO I = 1, IM
         XERR(I,J) = XERR(I,J) + BX(I,J) - BX(I,J+1)
      ENDDO
      ENDDO

      ! Return to calling program
      END SUBROUTINE POLFLT

!------------------------------------------------------------------------------

      SUBROUTINE DIAG_FLUX( IC, FX, FX1_TP, FY,  FY1_TP, 
     &                          FZ, FZ1_TP, NDT, ACOSP )    
!
!******************************************************************************
!  Subroutine DIAG_FLUX archives the mass fluxes in TPCORE version 7.1.
!  (bey, bmy, 9/20/00, 7/21/05)
!
!  Arguments as Input:
!  ============================================================================
!  (1  ) IC        (INTEGER) : Current tracer # 
!  (2,3) FX,FX1_TP (REAL*8 ) : Flux into the west side of grid box (I,J,K) 
!  (4,5) FY,FY1_TP (REAL*8 ) : Flux into the south side of grid box (I,J,K)
!  (6,7) FZ,FZ1_TP (REAL*8 ) : Flux into top of grid box (I,J,K) 
!  (8  ) NDT       (INTEGER) : Dynamic timestep in seconds
!  (9  ) ACOSP     (INTEGER) : Inverse cosine at latitude (J)
!
!  Included via header files:
!  ============================================================================
!  (1  ) DXYP   (REAL*8) : Surface area of grid box [m2]
!  (2  ) g0_100 (REAL*8) : The value 100 / 9.8
!
!  Diagnostics archived:
!  ============================================================================
!  (1  ) ND24 : Eastward flux of tracer in kg/s
!  (2  ) ND25 : Westward flux of tracer in kg/s
!  (3  ) ND26 : Upward   flux of tracer in kg/s
!
!  NOTES:
!  (1 ) Original code & algorithm is from Isabelle Bey, as installed in
!        TPCORE v. 4.1 (1998, 1999)
!  (2 ) DXYP is of dimension JJPAR, so reference it by DXYP(JREF), 
!        where JREF = J + J0. (bmy, 9/28/00)
!  (3 ) Add parallel processor directives to do-loops (bmy, 9/29/00)
!  (4 ) Archive CO budget array TCO for CO-OH run (bnd, bmy, 10/16/00)
!  (5 ) Also archive X-trop flux for CH4 simulation in TCH4 (bmy, 1/17/01)
!  (6 ) Added to "tpcore_mod.f" (bmy, 7/16/01)
!  (7 ) Now replace DXYP(JREF) with routine GET_AREA_M2 of "grid_mod.f".  
!        Also remove all references to JREF.  (bmy, 2/11/03)
!  (8 ) Now references TCVV and ITS_A_CH4_SIM from "tracer_mod.f" 
!        (bmy, 7/20/04)
!  (9 ) Remove references obsolete to CO-OH param code (bmy, 6/24/05)
!  (10) Bug fix: FX should be dimensioned with IIPAR+1 and FZ should be
!        dimensioned with LLPAR+1 (bmy, 7/21/05)
!  01 Mar 2012 - R. Yantosca - Now use GET_AREA_CM2(I,J,L) from grid_mod.F90
!******************************************************************************
!
      ! References to F90 modules
      USE DIAG_MOD,       ONLY : MASSFLEW, MASSFLNS, MASSFLUP
      USE GLOBAL_CH4_MOD, ONLY : XNUMOL_CH4, TCH4
      USE GRID_MOD,       ONLY : GET_AREA_M2
      USE TRACER_MOD,     ONLY : ITS_A_CH4_SIM, TCVV

      USE CMN_SIZE_MOD         ! Size parameters
      USE CMN_DIAG_MOD         ! Diagnostic switches
      USE CMN_GCTM_MOD         ! g0_100

      IMPLICIT NONE
#     include "define.h"

      ! Arguments
      INTEGER, INTENT(IN) :: IC, NDT
      REAL*8,  INTENT(IN) :: FX(IIPAR+1,JJPAR,LLPAR) 
      REAL*8,  INTENT(IN) :: FX1_TP(IIPAR,JJPAR,LLPAR)
      REAL*8,  INTENT(IN) :: FY(IIPAR,JJPAR,LLPAR)
      REAL*8,  INTENT(IN) :: FY1_TP(IIPAR,JJPAR,LLPAR)
      REAL*8,  INTENT(IN) :: FZ(IIPAR,JJPAR,LLPAR+1) 
      REAL*8,  INTENT(IN) :: FZ1_TP(IIPAR,JJPAR,LLPAR)
      REAL*8,  INTENT(IN) :: ACOSP(JJPAR)
       
      ! Local variables
      INTEGER             :: I, J, K, K2
      REAL*8              :: DTC, DTDYN, AREA_M2

      !=================================================================
      ! DIAG_FLUX begins here!
      ! 
      ! FX, FX1_TP, FY, FY1_TP, FZ, FZ1_TP have units of [mb/timestep].
      ! 
      ! To get tracer fluxes in kg/s :
      ! * (100./9.8)                 => kg/m2
      ! * DXYP(J)/(DTDYN * TCVV(IC)) => kg/s
      ! 
      ! Direction of the fluxes : 
      ! ----------------------------------------------------------------
      ! FX(I,J,K) => flux coming into the west edge of the box I 
      !              (from I-1 to I).
      !           => a positive flux goes from west to east.
      ! 
      ! FY(I,J,K) => flux coming into the south edge of the box J 
      !              (from J to J-1).
      !           => a positive flux goes from south to north     
      !              (from J-1 to J)
      ! 
      ! FZ(I,J,K) => flux coming down into the box k.
      !           => a positive flux goes down.
      !=================================================================

      ! DTDYN = double precision value for NDT, the dynamic timestep
      DTDYN = DBLE( NDT )

      !=================================================================
      ! ND24 Diagnostic: Eastward flux of tracer in [kg/s]
      !=================================================================
      IF ( ND24 > 0 ) THEN

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, K, K2, AREA_M2, DTC )
         DO K = 1, LLPAR

            ! K  is the vertical index down from the atmosphere top downwards
            ! K2 is the vertical index up from the surface
            K2 = LLPAR - K + 1

            DO J = 1, JJPAR

               DO I = 1, IIPAR
                  
                  ! Grid box surface area [m2]
                  AREA_M2 = GET_AREA_M2( I, J, K )

                  DTC = ( FX(I,J,K) + FX1_TP(I,J,K) ) * 
     &                  ( g0_100    * AREA_M2       ) /
     &                  ( TCVV(IC)  * DTDYN         )

                  MASSFLEW(I,J,K2,IC) = MASSFLEW(I,J,K2,IC) + DTC
               ENDDO
            ENDDO
         ENDDO
!$OMP END PARALLEL DO

      ENDIF

      !=================================================================
      ! ND25 Diagnostic: Northward flux of tracer in [kg/s]
      !=================================================================
      IF ( ND25 > 0 ) THEN
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, K, K2, AREA_M2, DTC )
         DO K = 1, LLPAR

            ! K  is the vertical index down from the atmosphere top downwards
            ! K2 is the vertical index up from the surface
            K2 = LLPAR - K + 1

            DO J = 1, JJPAR
               
               DO I = 1, IIPAR

                  ! Grid box surface area [m2]
                  AREA_M2 = GET_AREA_M2( I, J, K )

                  DTC = ( FY(I,J,K) + FY1_TP(I,J,K)    ) * 
     &                  ( ACOSP(J)  * g0_100 * AREA_M2 ) /
     &                  ( TCVV(IC)  * DTDYN            )

                  ! Contribution for CH4 run (bmy, 1/17/01)
                  IF ( ITS_A_CH4_SIM() ) THEN
                     TCH4(I,J,K,10) = TCH4(I,J,K,10) + 
     &                                ( DTC * DTDYN * XNUMOL_CH4 )
                  ENDIF

                  MASSFLNS(I,J,K2,IC) = MASSFLNS(I,J,K2,IC) + DTC
               ENDDO
            ENDDO
         ENDDO
!$OMP END PARALLEL DO

      ENDIF

      !=================================================================
      ! ND26 Diagnostic : Upward flux of tracer in [kg/s] 
      !=================================================================
      IF ( ND26 > 0 ) THEN
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, K, K2, DTC )
         DO K = 1, LLPAR

            ! K  is the vertical index down from the atmosphere top downwards
            ! K2 is the vertical index up from the surface
            K2 = LLPAR - K + 1

            DO J = 1, JJPAR

               DO I = 1, IIPAR

                  ! Grid box surface area [m2]
                  AREA_M2 = GET_AREA_M2( I, J, K )

                  DTC = ( FZ(I,J,K) + FZ1_TP(I,J,K) ) *
     &                  ( g0_100    * AREA_M2       ) /
     &                  ( TCVV(IC)  * DTDYN         )

                  MASSFLUP(I,J,K2,IC) = MASSFLUP(I,J,K2,IC) + DTC
               ENDDO
            ENDDO
         ENDDO
!$OMP END PARALLEL DO

      ENDIF

      ! Return to calling program
      END SUBROUTINE DIAG_FLUX

!------------------------------------------------------------------------------

      END MODULE TPCORE_MOD
