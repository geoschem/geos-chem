! $Id: pjc_pfix_geos5_window_mod.f,v 1.1 2009/09/16 14:06:15 bmy Exp $
      MODULE PJC_PFIX_GEOSFP_WINDOW_MOD
!
!******************************************************************************
!  Module PJC_PFIX_GEOS5_WINDOW_MOD contains routines which implements the 
!  Philip Cameron-Smith pressure fixer.  Specially modified for the GEOS-5
!  nested grid simulation. (yxw, dan, bmy, 11/6/08)
! 
!  Module Variables:
!  ============================================================================
!  (1 ) AI          (REAL*8 )  : Vertical coord "A" for hybrid grid [hPa]
!  (2 ) BI          (REAL*8 )  : Vertical coord "B" for hybrid grid [unitless]
!  (3 ) CLAT_FV     (REAL*8 )  : Grid box center latitude [radians] 
!  (4 ) COSE_FV     (REAL*8 )  : COSINE of grid box edge latitudes [radians]
!  (5 ) COSP_FV     (REAL*8 )  : COSINE of grid box center latitudes [radians]
!  (6 ) DAP         (REAL*8 )  : Delta-A vertical coordinate [hPa]
!  (7 ) DBK         (REAL*8 )  : Delta-B vertical coordinate [unitless]
!  (8 ) DLAT_FV     (REAL*8 )  : Latitude extent of grid boxes [radians]
!  (9 ) ELAT_FV     (REAL*8 )  : Grid box edge latitudes [radians]
!  (10) GEOFAC      (REAL*8 )  : Geometric factor for N-S advection
!  (11) GW_FV       (REAL*8 )  : Diff of SINE btw grid box lat edges [unitless]
!  (12) MCOR        (REAL*8 )  : Grid box surface areas [m2]
!  (13) REL_AREA    (REAL*8 )  : Relative surface area of grid box [fraction]
!  (14) RGW_FV      (REAL*8 )  : Reciprocal of GW_FV [radians
!  (15) SINE_FV     (REAL*8 )  : SINE of lat at grid box edges [unitless]
!  (16) GEOFAC_PC   (REAL*8 )  : Geometric factor for N-S advection @ poles
!  (17) DLON_FV     (REAL*8 )  : Longitude extent of a grid box [radians]
!  (18) LOC_PROC    (REAL*8 )  : Local processor number
!  (19) PR_DIAG     (LOGICAL)  : Flag for printing diagnostic message
!  (20) IMP_NBORDER (INTEGER)  : Used for ghost zones for MPI ??? 
!  (21) I1_GL       (INTEGER)  : ind of 1st  global lon       (no ghost zones)
!  (22) I2_GL       (INTEGER)  : ind of last global lon       (no ghost zones)
!  (23) JU1_GL      (INTEGER)  : ind of 1st  global "u" lat   (no ghost zones)
!  (24) JV1_GL      (INTEGER)  : ind of 1st  global "v" lat   (no ghost zones)
!  (25) J2_GL       (INTEGER)  : ind of last global "u&v" lat (no ghost zones)
!  (26) K1_GL       (INTEGER)  : ind of 1st  global alt       (no ghost zones)
!  (27) K2_GL       (INTEGER)  : ind of last global alt       (no ghost zones)
!  (28) ILO_GL      (INTEGER)  : I1_GL  - IMP_NBORDER        (has ghost zones)
!  (29) IHI_GL      (INTEGER)  : I2_GL  + IMP_NBORDER        (has ghost zones)
!  (30) JULO_GL     (INTEGER)  : JU1_GL - IMP_NBORDER        (has ghost zones)
!  (31) JVLO_GL     (INTEGER)  : JV1_GL - IMP_NBORDER        (has ghost zones)
!  (32) JHI_GL      (INTEGER)  : J2_GL  + IMP_NBORDER        (has ghost zones)
!  (33) I1          (INTEGER)  : ind of first local lon       (no ghost zones)
!  (34) I2          (INTEGER)  : ind of last  local lon       (no ghost zones)
!  (35) JU1         (INTEGER)  : ind of first local "u" lat   (no ghost zones)
!  (36) JV1         (INTEGER)  : ind of first local "v" lat   (no ghost zones)
!  (37) J2          (INTEGER)  : ind of last  local "u&v" lat (no ghost zones)
!  (38) K1          (INTEGER)  : index of first local alt     (no ghost zones)
!  (39) K2          (INTEGER)  : index of last  local alt     (no ghost zones)
!  (40) ILO         (INTEGER)  : I1  - IMP_NBORDER           (has ghost zones)
!  (41) IHI         (INTEGER)  : I2  + IMP_NBORDER           (has ghost zones)
!  (42) JULO        (INTEGER)  : JU1 - IMP_NBORDER           (has ghost zones)
!  (43) JVLO        (INTEGER)  : JV1 - IMP_NBORDER           (has ghost zones)
!  (44) JHI         (INTEGER)  : J2  + IMP_NBORDER           (has ghost zones)
!
!  Module Routines:
!  ============================================================================
!  (1 ) DO_PJC_PFIX            : Driver for Phil Cameron-Smith Pressure Fixer
!  (2 ) CHECK_TOTAL_MASS       : Prints total air mass and tracer masses
!  (3 ) CALC_PRESSURE          : Computes new pressures
!  (4 ) CALC_ADVECTION_FACTORS : Computes surface area factors for P-fixer
!  (5 ) ADJUST_PRESS           : Pressure-fixer routine from GMI/LLNL
!  (6 ) INIT_PRESS_FIX         : Pressure-fixer routine from GMI/LLNL
!  (7 ) DO_PRESS_FIX_LLNL      : Pressure-fixer routine from GMI/LLNL
!  (8 ) AVERAGE_PRESS_POLES    : Pressure-fixer routine from GMI/LLNL
!  (9 ) CONVERT_WINDS          : Pressure-fixer routine from GMI/LLNL
!  (10) CALC_HORIZ_MASS_FLUX   : Pressure-fixer routine from GMI/LLNL
!  (11) CALC_DIVERGENCE        : Pressure-fixer routine from GMI/LLNL
!  (12) SET_PRESS_TERMS        : Pressure-fixer routine from GMI/LLNL
!  (13) DO_DIVERGENCE_POLE_SUM : Pressure-fixer routine from GMI/LLNL
!  (14) XPAVG                  : Overwrites a 1-D vector w/ its avg value
!  (15) INIT_PJC_PFIX          : Initializes and allocates module variables
!  (16) CLEANUP_PJC_PFIX       : Deallocates all module variables
!
!  GEOS-CHEM modules referenced by tpcore_call_mod.f
!  ============================================================================
!  (1 ) error_mod.f    : Module containing I/O error/NaN check routines
!  (2 ) grid_mod.f     : Module containing horizontal grid information
!  (3 ) pressure_mod.f : Module containing routines to compute P(I,J,L)
!
!  NOTES:
!  (1 ) Adapted from "pjc_pfix_mod.f" (bmy, 11/6/08)
!  (2 ) Nested grids do not have polar cap and no periodicity at the edge.
!       Need to restrain the pressure to fixer to an inner window (-1 box in
!       all directions) because edges mass fluxes and divergence not known.
!       (lzh, ccc, 8/3/10)
!  01 Mar 2012 - R. Yantosca - Now reference new grid_mod.F90
!  20 Aug 2013 - R. Yantosca - Removed "define.h", this is now obsolete
!  26 Sep 2013 - R. Yantosca - Renamed to pjc_pfix_geosfp_window_mod.F
!******************************************************************************
!
      IMPLICIT NONE

      !=================================================================
      ! MODULE PRIVATE DECLARATIONS -- keep certain internal variables 
      ! and routines from being seen outside "pjc_pfix_mod.f"
      !=================================================================
      
      ! Make everything PRIVATE ...
      PRIVATE

      ! ... except these routines
      PUBLIC :: CLEANUP_PJC_PFIX_GEOSFP_WINDOW
      PUBLIC :: DO_PJC_PFIX_GEOSFP_WINDOW

      !=================================================================
      ! MODULE VARIABLES
      !=================================================================

      ! Allocatable arrays
      REAL*8, ALLOCATABLE :: AI(:)
      REAL*8, ALLOCATABLE :: BI(:)
      REAL*8, ALLOCATABLE :: CLAT_FV(:)
      REAL*8, ALLOCATABLE :: COSE_FV(:)
      REAL*8, ALLOCATABLE :: COSP_FV(:)
      REAL*8, ALLOCATABLE :: DAP(:)
      REAL*8, ALLOCATABLE :: DBK(:)
      REAL*8, ALLOCATABLE :: DLAT_FV(:)
      REAL*8, ALLOCATABLE :: ELAT_FV(:)
      REAL*8, ALLOCATABLE :: GEOFAC(:)
      REAL*8, ALLOCATABLE :: GW_FV(:)
      REAL*8, ALLOCATABLE :: MCOR(:,:)
      REAL*8, ALLOCATABLE :: REL_AREA(:,:)
      REAL*8, ALLOCATABLE :: RGW_FV(:)
      REAL*8, ALLOCATABLE :: SINE_FV(:)
 
      ! Scalar variables
      LOGICAL             :: PR_DIAG
      INTEGER             :: LOC_PROC
      REAL*8              :: GEOFAC_PC
      REAL*8              :: DLON_FV

      ! Dimensions for GMI code (from "imp_dims")
      INTEGER             :: IMP_NBORDER
      INTEGER             :: I1_GL,  I2_GL,   JU1_GL,  JV1_GL            
      INTEGER             :: J2_GL,  K1_GL,   K2_GL,   ILO_GL    
      INTEGER             :: IHI_GL, JULO_GL, JVLO_GL, JHI_GL    
      INTEGER             :: I1,     I2,      JU1,     JV1       
      INTEGER             :: J2,     K1,      K2,      ILO       
      INTEGER             :: IHI,    JULO,    JVLO,    JHI       
      INTEGER             :: ILAT,   ILONG,   IVERT,   J1P       
      INTEGER             :: J2P       

      ! Dimensions for nested grids 
      INTEGER             :: I1_W,     I2_W,      JU1_W       
      INTEGER             :: J2_W,     J1P_W,     J2P_W
      INTEGER             :: BUFF_SIZE

      !=================================================================
      ! MODULE ROUTINES -- follow below the "CONTAINS" statement
      !=================================================================
      CONTAINS

!------------------------------------------------------------------------------

      SUBROUTINE DO_PJC_PFIX_GEOSFP_WINDOW( D_DYN, P1,   P2, 
     &                                     UWND,  VWND, XMASS, YMASS )
!
!******************************************************************************
!  Subroutine DO_PJC_PFIX is the driver routine for the Philip Cameron-Smith
!  pressure fixer for the GEOS-4/fvDAS transport scheme. 
!  (bdf, bmy, 5/8/03, 3/5/07)
!
!  We assume that the winds are on the A-GRID, since this is the input that 
!  the GEOS-4/fvDAS transport scheme takes. (bdf, bmy, 5/8/03)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) D_DYN (REAL*8) : Dynamic timestep [s]
!  (2 ) P1    (REAL*8) : True PSurface at middle of dynamic timestep [hPa]
!  (3 ) P2    (REAL*8) : True PSurface at end    of dynamic timestep [hPa]
!  (4 ) UWND  (REAL*8) : Zonal (E-W) wind [m/s]
!  (5 ) VWND  (REAL*8) : Meridional (N-S) wind [m/s]
!
!  Arguments as Input:
!  ============================================================================
!  (6 ) XMASS (REAL*8) : E-W mass fluxes [kg/s ??]
!  (7 ) YMASS (REAL*8) : N-S mass fluxes [kg/s ??]
!  
!  NOTES:
!  (1 ) Now P1 and P2 are "true" surface pressures, and not PS-PTOP.  If using
!        this P-fixer w/ GEOS-3 winds, pass true surface pressure to this
!        routine. (bmy, 10/27/03)
!  (2 ) Now define P2_TMP array for passing to ADJUST_PRESS (yxw, bmy, 3/5/07)
!  20 Aug 2013 - R. Yantosca - Removed "define.h", this is now obsolete
!  26 Sep 2013 - R. Yantosca - Renamed to DO_PJC_PFIX_GEOSFP_WINDOW
!******************************************************************************
!
      USE CMN_SIZE_MOD    ! Size parameters
      USE CMN_GCTM_MOD    ! Physical constants

      IMPLICIT NONE

      ! Arguments
      REAL*8,  INTENT(IN)  :: D_DYN
      REAL*8,  INTENT(IN)  :: P1(IIPAR,JJPAR)
      REAL*8,  INTENT(IN)  :: P2(IIPAR,JJPAR)
      REAL*8,  INTENT(IN)  :: UWND(IIPAR,JJPAR,LLPAR)
      REAL*8,  INTENT(IN)  :: VWND(IIPAR,JJPAR,LLPAR)
      REAL*8,  INTENT(OUT) :: XMASS(IIPAR,JJPAR,LLPAR)
      REAL*8,  INTENT(OUT) :: YMASS(IIPAR,JJPAR,LLPAR)

      ! Local variables
      LOGICAL, SAVE        :: FIRST = .TRUE.
      INTEGER              :: I, J
      REAL*8               :: P2_TMP(IIPAR,JJPAR)
 
      ! Parameters
      LOGICAL, PARAMETER   :: INTERP_WINDS     = .TRUE.  ! winds are interp'd 
      INTEGER, PARAMETER   :: MET_GRID_TYPE    = 0       ! A-GRID
      INTEGER, PARAMETER   :: ADVEC_CONSRV_OPT = 0       ! 2=floating pressure 
      INTEGER, PARAMETER   :: PMET2_OPT        = 1       ! leave at 1 
      INTEGER, PARAMETER   :: PRESS_FIX_OPT    = 1       ! Turn on P-Fixer

      !=================================================================
      ! DO_PJC_PFIX begins here!
      !=================================================================

      ! Initialize on first call
      IF ( FIRST ) THEN

         ! Initialize/allocate module variables 
         CALL INIT_PJC_PFIX

         ! Calculate advection surface-area factors
         CALL CALC_ADVECTION_FACTORS( MCOR, REL_AREA, GEOFAC, GEOFAC_PC)

         ! Reset first-time flag
         FIRST = .FALSE.
      ENDIF

      ! Copy P2 into P2_TMP (yxw, bmy, 3/5/07)
      P2_TMP = P2
 
      ! Call PJC pressure fixer w/ the proper arguments
      ! NOTE: P1 and P2 are now "true" surface pressure, not PS-PTOP!!!
      CALL ADJUST_PRESS( 'GEOS-CHEM',        INTERP_WINDS,  
     &                   .TRUE.,             MET_GRID_TYPE, 
     &                   ADVEC_CONSRV_OPT,   PMET2_OPT, 
     &                   PRESS_FIX_OPT,      D_DYN, 
     &                   GEOFAC_PC,          GEOFAC, 
     &                   COSE_FV,            COSP_FV, 
     &                   REL_AREA,           DAP, 
     &                   DBK,                P1, 
     &                   P2_TMP,             P2_TMP,            
     &                   UWND,               VWND, 
     &                   XMASS,              YMASS )

      ! Return to calling program
      END SUBROUTINE DO_PJC_PFIX_GEOSFP_WINDOW

!------------------------------------------------------------------------------

      SUBROUTINE CALC_PRESSURE( XMASS, YMASS, RGW_FV, PS_NOW, PS_AFTER )
!
!******************************************************************************
!  Subroutine CALC_PRESSURE recalculates the new surface pressure from the
!  adjusted air masses XMASS and YMASS.  This is useful for debugging 
!  purposes. (bdf, bmy, 5/8/03)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) XMASS    (REAL*8) : E-W mass flux from pressure fixer
!  (2 ) YMASS    (REAL*8) : N-S mass flux from pressure fixer
!  (3 ) RGW_FV   (REAL*8) : 1 / ( SINE(J+1) - SINE(J) ) -- latitude factor
!  (4 ) PS_NOW   (REAL*8) : Surface pressure - PTOP at current time
!
!  Arguments as Output:
!  ============================================================================
!  (5 ) PS_AFTER (REAL*8) : Surface pressure - PTOP adjusted by P-fixer
!
!  NOTES:
!  20 Aug 2013 - R. Yantosca - Removed "define.h", this is now obsolete
!******************************************************************************
!
      USE CMN_SIZE_MOD  ! Size parameters
      USE CMN_MOD       ! NTRACE, LPRT, LWINDO

      IMPLICIT NONE


      ! Arguments
      REAL*8, INTENT(IN)  :: XMASS(IIPAR,JJPAR,LLPAR)
      REAL*8, INTENT(IN)  :: YMASS(IIPAR,JJPAR,LLPAR)
      REAL*8, INTENT(IN)  :: PS_NOW(IIPAR,JJPAR)
      REAL*8, INTENT(IN)  :: RGW_FV(JJPAR)
      REAL*8, INTENT(OUT) :: PS_AFTER(IIPAR,JJPAR)

      ! Local variables
      INTEGER             :: I, J, L
      REAL*8              :: DELP(IIPAR,JJPAR,LLPAR)
      REAL*8              :: DELP1(IIPAR,JJPAR,LLPAR)
      REAL*8              :: PE(IIPAR,LLPAR+1,JJPAR)

      !=================================================================
      ! CALC_PRESSURE begins here!
      !=================================================================
      DO L = 1, LLPAR
      DO J = 1, JJPAR
      DO I = 1, IIPAR
         DELP1(I,J,L) = DAP(L) + ( DBK(L) * PS_NOW(I,J) )
      ENDDO
      ENDDO
      ENDDO

      DO L = 1, LLPAR
         DO J = 2, JJPAR-1

            DO I =1, IIPAR-1
               DELP(I,J,L) = DELP1(I,J,L) + 
     &                       XMASS(I,J,L) - XMASS(I+1,J,L) + 
     &                     ( YMASS(I,J,L) - YMASS(I,J+1,L) ) * RGW_FV(J)
            ENDDO
            
            DELP(IIPAR,J,L) = 
     &          DELP1(IIPAR,J,L) +
     &          XMASS(IIPAR,J,L) - XMASS(1,J,L) +
     &        ( YMASS(IIPAR,J,L) - YMASS(IIPAR,J+1,L) ) * RGW_FV(J)
         ENDDO

         DO I = 1, IIPAR
            DELP(I,1,L) = DELP1(I,1,L) - YMASS(I,2,L) * RGW_FV(1)
         ENDDO

         ! Compute average
         CALL XPAVG( DELP(1,1,L), IIPAR )

         DO I = 1, IIPAR
            DELP(I,JJPAR,L) = DELP1(I,JJPAR,L) + 
     &                        YMASS(I,JJPAR,L) * RGW_FV(JJPAR)
         ENDDO

         ! Compute average
         CALL XPAVG( DELP(1,JJPAR,L), IIPAR )
      ENDDO

      !=================================================================
      ! Make the pressure
      !=================================================================
      DO J = 1, JJPAR
         DO I = 1, IIPAR
            PE(I,1,J) = PTOP
         ENDDO

         DO L = 1,LLPAR
            DO I = 1,IIPAR
               PE(I,L+1,J) = PE(I,L,J) + DELP(I,J,L)
            ENDDO
         ENDDO

         DO I = 1,IIPAR
            PS_AFTER(I,J) = PE(I,LLPAR+1,J)
         ENDDO
      ENDDO

      ! Return to calling program
      END SUBROUTINE CALC_PRESSURE

!------------------------------------------------------------------------------

      SUBROUTINE Calc_Advection_Factors
     &  (mcor, rel_area, geofac, geofac_pc)
!
!******************************************************************************
!
!  ROUTINE
!    Calc_Advection_Factors 
!  
!  AUTHORS
!    Philip Cameron-Smith and John Tannahill, GMI project @ LLNL (2003)
!  
!  DESCRIPTION 
!      This routine calculates the relative area of each grid
!      box, and the geometrical factors used by this modified version of
!      TPCORE. These geomoetrical DO assume that the space is regularly
!      grided, but do not assume any link between the surface area and
!      the linear dimensions.
!  
!  ARGUMENTS
!    mcor      : area of grid box (m^2)
!    rel_area  : relative surface area of grid box (fraction)
!    geofac    : geometrical factor for meridional advection; geofac uses
!                correct spherical geometry, and replaces acosp as the
!                meridional geometrical factor in tpcore
!    geofac_pc : special geometrical factor (geofac) for Polar cap
!  
!  NOTES:
!  (1 ) Now reference PI from "CMN_GCTM" for consistency.  Also force
!        double-precision with the "D" exponent. (bmy, 5/6/03)
!  (2 ) New definition for the geometric factor. (lzh, ccc, 8/3/10)
!******************************************************************************
!
      USE CMN_SIZE_MOD
      USE CMN_GCTM_MOD

      implicit none

      !----------------------
      !Argument declarations.
      !----------------------

      real*8, intent(IN)  :: mcor    (i1_gl :i2_gl, ju1_gl:j2_gl)
      real*8, intent(OUT) :: rel_area(i1_gl :i2_gl, ju1_gl:j2_gl)
      real*8, intent(OUT) :: geofac  (ju1_gl:j2_gl)
      real*8, intent(OUT) :: geofac_pc

      !----------------------
      !Variable declarations.
      !----------------------

      integer :: ij
 
! Variables not used. (ccc, 8/3/10)
!      real*8  :: dp           ! spacing in latitude (rad)
!      real*8  :: ri2_gl
!      real*8  :: rj2m1
      real*8  :: total_area

      !----------------
      !Begin execution.
      !----------------

! Not used. (ccc, 8/3/10)
!      ri2_gl = i2_gl

      !---------------------------------
      !Set the relative area (rel_area).
      !---------------------------------

      total_area = Sum (mcor(:,:))

      rel_area(:,:) = mcor(:,:) / total_area


      !---------------------------------------------------------
      !Calculate geometrical factor for meridional advection.
      !Note that it is assumed that all grid boxes in a latitude
      !band are the same.
      !---------------------------------------------------------

! Not used for nested grids. (ccc, 8/3/10)
!      rj2m1 = j2_gl - 1
!      dp    = PI / 360D0 

! The total area does not cover the full globe so use an other definition
! for the geometric factor. (lzh, ccc, 8/3/10)
      do ij = ju1_gl, j2_gl
!        geofac(ij) = dp / (2.0d0 * rel_area(1,ij) * ri2_gl)
        geofac(ij) = 1.d0 / COSP_FV(ij)
      end do

! geofac_pc used only for polar cap so no need. (ccc, 8/3/10)
!      geofac_pc =
!     &  dp / (2.0d0 * Sum (rel_area(1,ju1_gl:ju1_gl+1)) * ri2_gl)

      ! Return to calling program
      END SUBROUTINE Calc_Advection_Factors

!------------------------------------------------------------------------------

      SUBROUTINE Adjust_Press
     &  (metdata_name_org, do_timinterp_winds, new_met_rec,
     &   met_grid_type, advec_consrv_opt, pmet2_opt, press_fix_opt,
     &   tdt, geofac_pc, geofac, cose, cosp, rel_area, dap, dbk,
     &   pctm1, pctm2, pmet2, uu, vv, xmass, ymass)
!
!******************************************************************************
! 
!  ROUTINE
!    Adjust_Press
! 
!  AUTHORS
!    Philip Cameron-Smith and John Tannahill, GMI project @ LLNL (2003)
!
!  DESCRIPTION
!    This routine initializes and calls the pressure fixer code.
! 
!  ARGUMENTS
!    metdata_name_org   : first  part of metdata_name, e.g., "NCAR"
!    do_timinterp_winds : time interpolate wind fields?
!    new_met_rec        : new met record?
!    met_grid_type      : met grid type, A or C
!    advec_consrv_opt   : advection_conserve option
!    pmet2_opt          : pmet2 option
!    press_fix_opt      : pressure fixer option
!    tdt       : model time step (s)
!    geofac_pc : special geometrical factor (geofac) for Polar cap
!    geofac    : geometrical factor for meridional advection; geofac uses
!                correct spherical geometry, and replaces acosp as the
!                meridional geometrical factor in tpcore
!    cose      : cosine of grid box edges
!    cosp      : cosine of grid box centers
!    rel_area  : relative surface area of grid box (fraction)
!    dap       : pressure difference across layer from (ai * pt) term (mb)
!    dbk       : difference in bi across layer - the dSigma term
!    pctm1     : CTM surface pressure at t1     (mb)
!    pctm2     : CTM surface pressure at t1+tdt (mb)
!    pmet2     : surface pressure     at t1+tdt (mb)
!    uu        : wind velocity, x direction at t1+tdt/2 (m/s)
!    vv        : wind velocity, y direction at t1+tdt/2 (m/s)
!    xmass     : horizontal mass flux in E-W direction  (mb)
!    ymass     : horizontal mass flux in N-S direction  (mb)
! 
!  NOTES:
!  (1 ) Now declare METDATA_NAME_ORG as CHARACTER(LEN=*) (bmy, 6/25/03)
!******************************************************************************
!
      implicit none

      !----------------------
      !Argument declarations.
      !----------------------

      CHARACTER(LEN=*) :: metdata_name_org
      logical :: do_timinterp_winds
      logical :: new_met_rec
      integer :: met_grid_type
      integer :: advec_consrv_opt
      integer :: pmet2_opt
      integer :: press_fix_opt
      real*8  :: tdt
      real*8  :: geofac_pc
      real*8  :: geofac  (ju1_gl:j2_gl)
      real*8  :: cose    (ju1_gl:j2_gl)
      real*8  :: cosp    (ju1_gl:j2_gl)
      real*8  :: dap     (k1:k2)
      real*8  :: dbk     (k1:k2)

      !-----------------------------------------------------------------
      !rel_area : relative surface area of grid box (fraction)
      !-----------------------------------------------------------------
 
      real*8 :: rel_area( i1_gl:i2_gl,   ju1_gl:j2_gl)
            
      !-----------------------------------------------------------------
      !pmet1  : Metfield surface pressure at t1 (mb)
      !pmet2  : Metfield surface pressure at t1+tdt (mb)
      !pctm1  : CTM surface pressure at t1 (mb)
      !pctm2  : CTM surface pressure at t1+tdt (mb)
      !-----------------------------------------------------------------

      REAL*8 :: 
     &   pmet2(ilo_gl:ihi_gl, julo_gl:jhi_gl),
     &   pctm1(ilo_gl:ihi_gl, julo_gl:jhi_gl),
     &   pctm2(ilo_gl:ihi_gl, julo_gl:jhi_gl)

      !-------------------------------------------------
      !uu : wind velocity, x direction at t1+tdt/2 (m/s)
      !vv : wind velocity, y direction at t1+tdt/2 (m/s)
      !-------------------------------------------------

      real*8 :: uu(ilo_gl:ihi_gl, julo_gl:jhi_gl, k1_gl:k2_gl)
      real*8 :: vv(ilo_gl:ihi_gl, julo_gl:jhi_gl, k1_gl:k2_gl)

      !--------------------------------------------------
      !xmass : horizontal mass flux in E-W direction (mb)
      !ymass : horizontal mass flux in N-S direction (mb)
      !--------------------------------------------------

      real*8 :: xmass(ilo_gl:ihi_gl, julo_gl:jhi_gl, k1_gl:k2_gl)
      real*8 :: ymass(ilo_gl:ihi_gl, julo_gl:jhi_gl, k1_gl:k2_gl)

      logical, save :: DO_ADJUST_PRESS_DIAG = .TRUE.


      !----------------------
      !Variable declarations.
      !----------------------

      logical, save :: first = .true.

      !--------------------------------------------------
      !dgpress   : global-pressure discrepancy
      !press_dev : RMS difference between pmet2 and pctm2
      !            (weighted by relative area)
      !--------------------------------------------------

      real*8  :: dgpress
      real*8  :: press_dev

      !-------------------------------------------------------------
      !dps : change of surface pressure from met field pressure (mb)
      !-------------------------------------------------------------

      real*8  :: dps(i1_gl:i2_gl, ju1_gl:j2_gl)

      !--------------------------------------------
      !dps_ctm : CTM surface pressure tendency (mb)
      !--------------------------------------------

      real*8 :: dps_ctm(i1_gl:i2_gl, ju1_gl:j2_gl)

      !---------------------------------------------------------------------
      !xmass_fixed : horizontal mass flux in E-W direction after fixing (mb)
      !ymass_fixed : horizontal mass flux in N-S direction after fixing (mb)
      !---------------------------------------------------------------------

      real*8  :: xmass_fixed(ilo_gl:ihi_gl, julo_gl:jhi_gl, k1:k2)
      real*8  :: ymass_fixed(ilo_gl:ihi_gl, julo_gl:jhi_gl, k1:k2)

      !-------------
      !Dummy indexes
      !-------------

      !integer :: ij, il

      !----------------
      !Begin execution.
      !----------------

      if (pr_diag) then
        Write (6, *) 'Adjust_Press called by ', loc_proc
      end if

      dps_ctm(:,:) = 0.0d0

      dgpress =  Sum ( (pmet2(i1_gl:i2_gl, ju1_gl:j2_gl) -
     &                  pctm1(i1_gl:i2_gl, ju1_gl:j2_gl)   )
     &             * rel_area(i1_gl:i2_gl, ju1_gl:j2_gl)     )

      if (pmet2_opt == 1) then
        pmet2(:,:) = pmet2(:,:) - dgpress
      end if

!### Debug
!###      if (DO_ADJUST_PRESS_DIAG) then
!###        Write (6, *) 'Global mean surface pressure change (mb) = ', 
!###     &                dgpress
!###      end if

       !===================
        call Init_Press_Fix
       !===================
     &    (metdata_name_org, met_grid_type, tdt, geofac_pc, geofac,
     &     cose, cosp, dap, dbk, dps, dps_ctm, rel_area, pctm1, pmet2,
     &     uu, vv, xmass, ymass)

        if (press_fix_opt == 1) then

         !======================
          call Do_Press_Fix_Llnl
         !======================
     &      (geofac_pc, geofac, dbk, dps, dps_ctm, rel_area,
     &       xmass, ymass, xmass_fixed, ymass_fixed )
          
          xmass(:,:,:) = xmass_fixed(:,:,:)
          ymass(:,:,:) = ymass_fixed(:,:,:)

        end if

        if ((advec_consrv_opt == 0) .or.
     &      (advec_consrv_opt == 1)) then

          dps_ctm(i1_gl:i2_gl, ju1_gl:j2_gl) =
     &      pmet2(i1_gl:i2_gl, ju1_gl:j2_gl) - 
     &      pctm1(i1_gl:i2_gl, ju1_gl:j2_gl) 

        !-----------------------------------------------
        !else if (advec_consrv_opt == 2) then do nothing
        !-----------------------------------------------

        end if


      pctm2(i1_gl:i2_gl, ju1_gl:j2_gl) =
     &      pctm1(i1_gl:i2_gl, ju1_gl:j2_gl) + 
     &    dps_ctm(i1_gl:i2_gl, ju1_gl:j2_gl)


      if (DO_ADJUST_PRESS_DIAG) then

        !-------------------------------------------------------
        !Calculate the RMS pressure deviation (diagnostic only).
        !-------------------------------------------------------

        press_dev = 
     &    Sqrt (Sum (((pmet2(i1_gl:i2_gl,ju1_gl:j2_gl) -
     &                 pctm2(i1_gl:i2_gl,ju1_gl:j2_gl))**2 *
     &                rel_area(i1_gl:i2_gl,ju1_gl:j2_gl))))

!### Debug
!###        Write (6, *) 'RMS deviation between pmet2 & pctm2 (mb) = ',
!###     &               press_dev

      end if

      ! Return to calling program
      END SUBROUTINE Adjust_Press

!------------------------------------------------------------------------------

      SUBROUTINE Init_Press_Fix
     &  (metdata_name_org, met_grid_type, tdt, geofac_pc, geofac,
     &   cose, cosp, dap, dbk, dps, dps_ctm, rel_area, pctm1, pmet2,
     &   uu, vv, xmass, ymass)
!
!******************************************************************************
! 
!  ROUTINE
!    Init_Press_Fix
!
!  AUTHORS
!    Philip Cameron-Smith and John Tannahill, GMI project @ LLNL (2003)
! 
!  DESCRIPTION
!    This routine initializes the pressure fixer.
! 
!  ARGUMENTS
!    metdata_name_org : first  part of metdata_name, e.g., "NCAR"
!    met_grid_type    : met grid type, A or C
!    tdt       : model time step (s)
!    geofac_pc : special geometrical factor (geofac) for Polar cap
!    geofac    : geometrical factor for meridional advection; geofac uses
!                correct spherical geometry, and replaces acosp as the
!                meridional geometrical factor in tpcore
!    cose      : cosine of grid box edges
!    cosp      : cosine of grid box centers
!    dap       : pressure difference across layer from (ai * pt) term (mb)
!    dbk       : difference in bi across layer - the dSigma term
!    dps       : change of surface pressure from met field pressure (mb)
!    dps_ctm   : sum over vertical of dpi calculated from original mass
!                fluxes (mb)
!    rel_area  : relative surface area of grid box (fraction)
!    pctm1     : CTM       surface pressure at t1      (mb)
!    pmet2     : met field surface pressure at t1+tdt  (mb)
!    uu        : wind velocity in E-W direction        (m/s)
!    vv        : wind velocity in N-S direction        (m/s)
!    xmass     : horizontal mass flux in E-W direction (mb)
!    ymass     : horizontal mass flux in N-S direction (mb)
!
!  NOTES: 
!  (1 ) Now declare METDATA_NAME_ORG as CHARACTER(LEN=*)
!******************************************************************************
! 
      implicit none


      !----------------------
      !Argument declarations.
      !----------------------

      CHARACTER(LEN=*) :: metdata_name_org
      integer :: met_grid_type
      real*8  :: geofac_pc
      real*8  :: tdt
      real*8  :: cose    (ju1_gl:j2_gl)
      real*8  :: cosp    (ju1_gl:j2_gl)
      real*8  :: geofac  (ju1_gl:j2_gl)
      real*8  :: dap     (k1:k2)
      real*8  :: dbk     (k1:k2)

      !-------------------------------------------------------------
      !dps : change of surface pressure from met field pressure (mb)
      !-------------------------------------------------------------

      real*8  :: dps(i1_gl:i2_gl, ju1_gl:j2_gl)

      !--------------------------------------------
      !dps_ctm : CTM surface pressure tendency (mb)
      !--------------------------------------------

      real*8 :: dps_ctm(i1_gl:i2_gl, ju1_gl:j2_gl)

      !-----------------------------------------------------------------
      !rel_area : relative surface area of grid box (fraction)
      !-----------------------------------------------------------------

      real*8 :: rel_area( i1_gl:i2_gl,   ju1_gl:j2_gl)
            
      !-----------------------------------------------------------------
      !pmet1  : Metfield surface pressure at t1 (mb)
      !pmet2  : Metfield surface pressure at t1+tdt (mb)
      !pctm1  : CTM surface pressure at t1 (mb)
      !pctm2  : CTM surface pressure at t1+tdt (mb)
      !-----------------------------------------------------------------

      REAL*8 :: 
     &   pmet2(ilo_gl:ihi_gl, julo_gl:jhi_gl),
     &   pctm1(ilo_gl:ihi_gl, julo_gl:jhi_gl),
     &   pctm2(ilo_gl:ihi_gl, julo_gl:jhi_gl)

      !--------------------------------------------------
      ! uu : wind velocity, x direction at t1+tdt/2 (m/s)
      ! vv : wind velocity, y direction at t1+tdt/2 (m/s)
      !--------------------------------------------------

      real*8 :: uu(ilo_gl:ihi_gl, julo_gl:jhi_gl, k1_gl:k2_gl)
      real*8 :: vv(ilo_gl:ihi_gl, julo_gl:jhi_gl, k1_gl:k2_gl)

      !--------------------------------------------------
      !xmass : horizontal mass flux in E-W direction (mb)
      !ymass : horizontal mass flux in N-S direction (mb)
      !--------------------------------------------------

      real*8 :: xmass(ilo_gl:ihi_gl, julo_gl:jhi_gl, k1_gl:k2_gl)
      real*8 :: ymass(ilo_gl:ihi_gl, julo_gl:jhi_gl, k1_gl:k2_gl)

      !----------------------
      !Variable declarations.
      !----------------------

      !--------------------------------------------------------------
      !dpi   : divergence at a grid point; used to calculate vertical
      !        motion (mb)
      !--------------------------------------------------------------

      real*8  :: dpi(i1:i2, ju1:j2, k1:k2)

      !---------------------------------------------------------------------
      !crx   : Courant number in E-W direction
      !cry   : Courant number in N-S direction
      !delp1 : pressure thickness, the psudo-density in a hydrostatic system
      !        at t1 (mb)
      !delpm : pressure thickness, the psudo-density in a hydrostatic system
      !        at t1+tdt/2 (approximate) (mb)
      !pu    : pressure at edges in "u"  (mb)
      !---------------------------------------------------------------------

      real*8  :: crx  (ilo:ihi, julo:jhi, k1:k2)
      real*8  :: cry  (ilo:ihi, julo:jhi, k1:k2)
      real*8  :: delp1(ilo:ihi, julo:jhi, k1:k2)
      real*8  :: delpm(ilo:ihi, julo:jhi, k1:k2)
      real*8  :: pu   (ilo:ihi, julo:jhi, k1:k2)

      !----------------
      !Begin execution.
      !----------------

      if (pr_diag) then
        Write (6,*) 'Init_Press_Fix called by ', loc_proc
      end if

! not treat poles (lzh, 07/20/2010)
!      !========================
!      call Average_Press_Poles
!      !========================
!     &  (rel_area, pctm1)
!
!      !========================
!      call Average_Press_Poles
!      !========================
!     &  (rel_area, pmet2)

      !-------------------------------------------------------------------
      !We need to calculate pressures at t1+tdt/2.  One ought to use pctm2
      !in the call to Set_Press_Terms, but since we don't know it yet, we
      !are forced to use pmet2.  This should be good enough because it is
      !only used to convert the winds to the mass fluxes, which is done
      !crudely anyway and the mass fluxes will still get fixed OK.
      !-------------------------------------------------------------------

      dps(i1:i2,ju1:j2) = pmet2(i1:i2,ju1:j2) - pctm1(i1:i2,ju1:j2)

      !====================
      call Set_Press_Terms
      !====================
     &  (dap, dbk, pctm1, pmet2, delp1, delpm, pu)


        !===================
        call Convert_Winds
        !===================
     &    (met_grid_type, tdt, cosp, crx, cry, uu, vv)
        

        !=========================
        call Calc_Horiz_Mass_Flux
        !=========================
     &    (cose, delpm, uu, vv, xmass, ymass, tdt, cosp)

      !====================
      call Calc_Divergence
      !====================
     &  (.false., geofac_pc, geofac, dpi, xmass, ymass)

        
      dps_ctm(i1:i2,ju1:j2) = Sum (dpi(i1:i2,ju1:j2,:), dim=3)

      ! Return to calling program
      END SUBROUTINE Init_Press_Fix

!------------------------------------------------------------------------------

      SUBROUTINE Do_Press_Fix_Llnl
     &  (geofac_pc, geofac, dbk, dps, dps_ctm, rel_area,
     &   xmass, ymass, xmass_fixed, ymass_fixed)
!
!******************************************************************************
! 
!  ROUTINE
!    Do_Press_Fix_Llnl
! 
!  AUTHORS
!    Philip Cameron-Smith and John Tannahill, GMI project @ LLNL (2003)
!
!  DESCRIPTION
!    This routine fixes the mass fluxes to match the met field pressure
!    tendency.
! 
!  ARGUMENTS
!    geofac_pc   : special geometrical factor (geofac) for Polar cap
!    geofac      : geometrical factor for meridional advection; geofac uses
!                  correct spherical geometry, and replaces acosp as the
!                  meridional geometrical factor in tpcore
!    dbk         : difference in bi across layer - the dSigma term
!    dps         : change of surface pressure from met field pressure (mb)
!    dps_ctm     : sum over vertical of dpi calculated from original mass
!                  fluxes (mb)
!    rel_area    : relative surface area of grid box (fraction)
!    xmass       : horizontal mass flux in E-W direction (mb)
!    ymass       : horizontal mass flux in N-S direction (mb)
!    xmass_fixed : horizontal mass flux in E-W direction after fixing (mb)
!    ymass_fixed : horizontal mass flux in N-S direction after fixing (mb)
! 
!  NOTES:
! 
!******************************************************************************
!
      implicit none

      !----------------------
      !Argument declarations.
      !----------------------

      real*8  :: geofac_pc
      real*8  :: geofac  (ju1_gl:j2_gl)
      real*8  :: dbk     (k1:k2)
      real*8  :: dps     (i1:i2, ju1:j2)
      real*8  :: dps_ctm (i1:i2, ju1:j2)
      real*8  :: rel_area(i1:i2, ju1:j2)
      real*8  :: xmass      (ilo:ihi, julo:jhi, k1:k2)
      real*8  :: ymass      (ilo:ihi, julo:jhi, k1:k2)
      real*8  :: xmass_fixed(ilo:ihi, julo:jhi, k1:k2)
      real*8  :: ymass_fixed(ilo:ihi, julo:jhi, k1:k2)

      !----------------------
      !Variable declarations.
      !----------------------

      integer :: il, ij, ik

      real*8  :: dgpress
      real*8  :: fxmean
      real*8  :: ri2

      real*8  :: fxintegral(i1:i2+1)

      real*8  :: mmfd(ju1:j2)
      real*8  :: mmf (ju1:j2)

      real*8  :: ddps(i1:i2, ju1:j2)

      !------------------------------------------------------------------------
      !dpi : divergence at a grid point; used to calculate vertical motion (mb)
      !------------------------------------------------------------------------

      real*8  :: dpi(i1:i2, ju1:j2, k1:k2)

      real*8  :: xcolmass_fix(ilo:ihi, julo:jhi)

      real*8  :: xx

      !----------------
      !Begin execution.
      !----------------

      if (pr_diag) then
        Write (6,*) 'Do_Press_Fix_Llnl called by ', loc_proc
      end if


      ri2 = i2_gl - 2 * BUFF_SIZE

      mmfd(:) = 0.0d0
      mmf(:) = 0d0

      xcolmass_fix(:,:)   = 0.0d0

      xmass_fixed (:,:,:) = xmass(:,:,:)
      ymass_fixed (:,:,:) = ymass(:,:,:)


      !------------------------------------------------------------
      !Calculate difference between GCM and LR predicted pressures.
      !------------------------------------------------------------

      ddps(:,:) = dps(:,:) - dps_ctm(:,:)


c     --------------------------------------
c     Calculate global-pressure discrepancy.
c     --------------------------------------

!      dgpress =
!     &  Sum (ddps(i1:i2,ju1:j2) * rel_area(i1:i2,ju1:j2))

      xx = sum(rel_area(i1_w:i2_w,j1p_w:j2p_w))
      dgpress =
     &  Sum (ddps(i1_w:i2_w,j1p_w:j2p_w) * 
     &  rel_area(i1_w:i2_w,j1p_w:j2p_w))
     &  / xx


      !----------------------------------------------------------
      !Calculate mean meridional flux divergence (df/dy).
      !Note that mmfd is actually the zonal mean pressure change,
      !which is related to df/dy by geometrical factors.
      !----------------------------------------------------------

      !------------------------
      !Handle non-Pole regions.
      !------------------------

! Work on the inner window only (lzh, ccc, 8/3/10)
!      do ij = j1p, j2p
!        mmfd(ij) = -(sum(ddps(:,ij)) / ri2 - dgpress) 
!      end do

      do ij = j1p_w, j2p_w
        mmfd(ij) = -(sum(ddps(i1_w:i2_w,ij)) / ri2 - dgpress) 
      end do

! No special case for poles, no poles. (ccc, 8/3/10)
!      !---------------------------------------------
!      !Handle poles.
!      !Note that polar boxes have all been averaged.
!      !---------------------------------------------
!
!       mmfd(ju1)   = -(ddps(1,ju1)   - dgpress)
!       mmfd(ju1+1) = -(ddps(1,ju1+1) - dgpress)
!       mmfd(j2-1)  = -(ddps(1,j2-1)  - dgpress)
!       mmfd(j2)    = -(ddps(1,j2)    - dgpress)


      !---------------------------------------------
      !Calculate mean meridional fluxes (cos(e)*fy).
      !---------------------------------------------

! Use geofac, no polar cap. (ccc, 8/3/10)
!       mmf(j1p) = mmfd(ju1) / geofac_pc
       mmf(j1p_w) = mmfd(ju1_w) / geofac(j1p_w)

! Work on inner domain. (ccc, 8/3/10)
!       do ij = j1p, j2p
       do ij = j1p_w, j2p_w-1
          mmf(ij+1) = mmf(ij) + mmfd(ij) / geofac(ij)
       end do


      !------------------------------------------------------------
      !Fix latitude bands.
      !Note that we don't need to worry about geometry here because
      !all boxes in a latitude band are identical.
      !Note also that fxintegral(i2+1) should equal fxintegral(i1),
      !i.e., zero.
      !------------------------------------------------------------

! Work on inner domain (ccc, 8/3/10)
!      do ij = j1p, j2p
      do ij = j1p_w, j2p_w

        fxintegral(:) = 0.0d0

!        do il = i1, i2
        do il = i1_w, i2_w
          fxintegral(il+1) =
     &      fxintegral(il) -
     &      (ddps(il,ij) - dgpress) -
     &      mmfd(ij)
        end do

        fxmean = Sum (fxintegral(i1+1:i2+1)) / ri2
!        fxmean = Sum (fxintegral(i1_w+1:i2_w+1)) / ri2

!        do il = i1, i2
        do il = i1_w, i2_w
          xcolmass_fix(il,ij) = fxintegral(il) - fxmean
        end do

      end do

      !-------------------------------------
      !Distribute colmass_fix's in vertical.
      !-------------------------------------

      do ik = k1, k2      
!        do ij = j1p, j2p
!          do il = i1, i2
        do ij = j1p_w, j2p_w
          do il = i1_w, i2_w

            xmass_fixed(il,ij,ik) = xmass(il,ij,ik) + 
     &                              xcolmass_fix(il,ij) * dbk(ik)

          end do
        end do
      end do

! Grid stops at j2p if nested domain (ccc, 8/3/10)
!      do ik = k1, k2      
!        do ij = j1p, j2p+1
!          do il = i1, i2
!
!            ymass_fixed(il,ij,ik) = ymass(il,ij,ik) +
!     &                              mmf(ij) * dbk(ik)
!
!          end do
!        end do
!      end do
      do ik = k1, k2      
!        do ij = j1p, j2p
!          do il = i1, i2
        do ij = j1p_w, j2p_w
          do il = i1_w, i2_w

            ymass_fixed(il,ij,ik) = ymass(il,ij,ik) +
     &                              mmf(ij) * dbk(ik)

          end do
        end do
      end do
 
      !====================
      call Calc_Divergence
      !====================
     &  (.false., geofac_pc, geofac, dpi, xmass_fixed, ymass_fixed)


      dps_ctm(i1:i2,ju1:j2) = Sum (dpi(i1:i2,ju1:j2,:), dim=3)

      ! Return to calling program
      END SUBROUTINE Do_Press_Fix_Llnl

!------------------------------------------------------------------------------

      SUBROUTINE Convert_Winds
     &  (igd, tdt, cosp, crx, cry, uu, vv)
!
!******************************************************************************
!
!  ROUTINE
!    Convert_Winds
!  
!  AUTHORS
!    Philip Cameron-Smith and John Tannahill, GMI project @ LLNL (2003)
!  
!  DESCRIPTION
!    This routine converts winds on A or C grid to Courant # on C grid.
!  
!  ARGUMENTS
!    igd  : A or C grid
!    tdt  : model time step (s)
!    cosp : cosine of grid box centers
!    crx  : Courant number in E-W direction
!    cry  : Courant number in N-S direction
!    uu   : wind velocity  in E-W direction at t1+tdt/2 (m/s)
!    vv   : wind velocity  in N-S direction at t1+tdt/2 (m/s)
!  
!  NOTES:
!  (1 ) Use GEOS-CHEM physical constants Re, PI to be consistent with other
!        usage everywhere (bmy, 5/5/03)
!
!******************************************************************************
!
      USE CMN_SIZE_MOD ! Size parameters
      USE CMN_GCTM_MOD ! Re, PI
 
      implicit none

      !----------------------
      !Argument declarations.
      !----------------------

      integer :: igd
      real*8  :: tdt
      real*8  :: cosp(ju1_gl:j2_gl)
      real*8  :: crx (ilo:ihi, julo:jhi, k1:k2)
      real*8  :: cry (ilo:ihi, julo:jhi, k1:k2)
      real*8  :: uu  (ilo:ihi, julo:jhi, k1:k2)
      real*8  :: vv  (ilo:ihi, julo:jhi, k1:k2)

      !----------------------
      !Variable declarations.
      !----------------------

      logical, save :: first = .true.

      integer :: il, ij

      !-------------------------------
      !dl : spacing in longitude (rad)
      !dp : spacing in latitude  (rad)
      !-------------------------------

      real*8  :: dl
      real*8  :: dp

      real*8  :: ri2
      real*8  :: rj2m1

      !------------------------
      !dtdy  : dt/dy      (s/m)
      !dtdy5 : 0.5 * dtdy (s/m)
      !------------------------

      real*8, save :: dtdy
      real*8, save :: dtdy5

      !------------------------
      !dtdx  : dt/dx      (s/m)
      !dtdx5 : 0.5 * dtdx (s/m)
      !------------------------

      real*8, allocatable, save :: dtdx (:)                   
      real*8, allocatable, save :: dtdx5(:)             

      !----------------
      !Begin execution.
      !----------------

      if (pr_diag) then
        Write (6, *) 'Convert_Winds called by ', loc_proc
      end if


      !==========
      if (first) then
      !==========

        first = .false.

        Allocate (dtdx (ju1_gl:j2_gl))
        Allocate (dtdx5(ju1_gl:j2_gl))

        dtdx = 0.0d0; dtdx5 = 0.0d0

        ri2   = i2_gl
        rj2m1 = j2_gl - 1

        dl    = 2.0d0 * PI / 1152D0  !(dan)
        dp    = PI /720D0            !(dan)

        dtdy  = tdt / (Re * dp)
        dtdy5 = 0.5d0 * dtdy

!-----lzh----------
!        dtdx (ju1_gl) = 0.0d0
!        dtdx5(ju1_gl) = 0.0d0
!
!        do ij = ju1_gl + 1, j2_gl - 1
!
!          dtdx (ij) = tdt / (dl * Re * cosp(ij))
!          dtdx5(ij) = 0.5d0 * dtdx(ij)
!
!        end do
!
!        dtdx (j2_gl)  = 0.0d0
!        dtdx5(j2_gl)  = 0.0d0

!-----------------------------------------------
! for nested NA or EA (lzh, 07/20/2010)
        do ij = ju1_gl, j2_gl
          dtdx (ij) = tdt / (dl * Re * cosp(ij))
          dtdx5(ij) = 0.5d0 * dtdx(ij)
        end do
!-----------------------------------------------

      end if


      !=============
      if (igd == 0) then  ! A grid.
      !=============

        do ij = ju1+1, j2-1
          do il = i1+1, i2
            crx(il,ij,:) =
     &        dtdx5(ij) *
     &        (uu(il,ij,:) + uu(il-1,ij,  :))
          end do
! No periodicity (ccc, 8/3/10)
!            crx(1,ij,:) =
!     &        dtdx5(ij) *
!     &        (uu(1,ij,:) + uu(i2,ij,  :))
        end do

        do ij = ju1+1, j2
          do il = i1, i2
            cry(il,ij,:) =
     &        dtdy5 *
     &        (vv(il,ij,:) + vv(il,  ij-1,:))
          end do
        end do


      !====
      else  ! C grid.
      !====

! No ghost zones. (ccc, 8/3/10)
!        do ij = ju1, j2
!          do il = i1, i2
        do ij = ju1+1, j2
          do il = i1+1, i2

            crx(il,ij,:) =
     &        dtdx(ij) * uu(il-1,ij,  :)

            cry(il,ij,:) =
     &        dtdy     * vv(il,  ij-1,:)

          end do
        end do

      end if

      ! Return to calling program
      END SUBROUTINE Convert_Winds

!------------------------------------------------------------------------------

      SUBROUTINE Calc_Horiz_Mass_Flux
     &  (cose, delpm, uu, vv, xmass, ymass, tdt, cosp)
!
!******************************************************************************
! 
!  ROUTINE
!    Calc_Horiz_Mass_Flux
! 
!  AUTHORS
!    Philip Cameron-Smith and John Tannahill, GMI project @ LLNL (2003)
! 
!  DESCRIPTION
!    This routine calculates the horizontal mass flux for non-GISS met data.
! 
!  ARGUMENTS
!    cose  : cosine of grid box edges
!    delpm : pressure thickness, the psudo-density in a hydrostatic system
!            at t1+tdt/2 (approximate) (mb)
!    crx   : Courant number in E-W direction
!    cry   : Courant number in N-S direction
!    pu    : pressure at edges in "u"  (mb)
!    xmass : horizontal mass flux in E-W direction (mb)
!    ymass : horizontal mass flux in N-S direction (mb)
! 
!  NOTES:
!  
!******************************************************************************
!
      USE CMN_SIZE_MOD ! Size parameters
      USE CMN_GCTM_MOD ! Re, Pi

      implicit none

      !----------------------
      !Argument declarations.
      !----------------------

      real*8  :: tdt
      real*8  :: cose (ju1_gl:j2_gl)
      real*8  :: cosp (ju1_gl:j2_gl)
      real*8  :: delpm(ilo:ihi, julo:jhi, k1:k2)
      real*8  :: uu  (ilo:ihi, julo:jhi, k1:k2)
      real*8  :: vv  (ilo:ihi, julo:jhi, k1:k2)
      real*8  :: xmass(ilo:ihi, julo:jhi, k1:k2)
      real*8  :: ymass(ilo:ihi, julo:jhi, k1:k2)

      !----------------------
      !Variable declarations.
      !----------------------

      integer :: ij
      integer :: il
      integer :: jst, jend
      real*8  :: dl
      real*8  :: dp

      real*8  :: ri2
      real*8  :: rj2m1
      real*8  :: factx
      real*8  :: facty

      !----------------
      !Begin execution.
      !----------------

      if (pr_diag) then
        Write (6,*) 'Calc_Horiz_Mass_Flux called by ', loc_proc
      end if

        ri2   = i2_gl
        rj2m1 = j2_gl - 1

        dl    = 2.0d0 * PI /1152D0   !(dan) 
        dp    = PI /720D0            !(dan) 

        facty  = 0.5d0 * tdt / (Re * dp)

      !-----------------------------------
      !Calculate E-W horizontal mass flux.
      !-----------------------------------

      do ij = ju1, j2

       factx = 0.5d0 * tdt / (dl * Re * cosp(ij))

       do il = i1+1, i2 
        xmass(il,ij,:) = factx *
     &    (uu(il,ij,:) * delpm(il,ij,:)+
     &     uu(il-1,ij,:) * delpm(il-1,ij,:))
       end do

! No periodicity. (ccc, 8/3/10)
!        xmass(i1,ij,:) = factx *
!     &    (uu(i1,ij,:) * delpm(i1,ij,:)+
!     &     uu(i2,ij,:) * delpm(i2,ij,:))

      end do


      !-----------------------------------
      !Calculate N-S horizontal mass flux.
      !-----------------------------------

      do ij = ju1+1, j2

         ymass(i1:i2,ij,:) = facty *
     &    cose(ij) * (vv(i1:i2,ij,:)*delpm(i1:i2,ij,:)+
     &    vv(i1:i2,ij-1,:)*delpm(i1:i2,ij-1,:))

      end do

      ! Return to calling program
      END SUBROUTINE Calc_Horiz_Mass_Flux

!------------------------------------------------------------------------------

      SUBROUTINE Calc_Divergence
     &  (do_reduction, geofac_pc, geofac, dpi, xmass, ymass)
!
!******************************************************************************
!
!  ROUTINE
!    Calc_Divergence
!
!  AUTHORS
!    Philip Cameron-Smith and John Tannahill, GMI project @ LLNL (2003)
!
!  DESCRIPTION
!    This routine calculates the divergence.
!
!  ARGUMENTS
!    do_reduction : set to false if called on Master;
!                   set to true  if called by Slaves
!    geofac_pc    : special geometrical factor (geofac) for Polar cap
!    geofac       : geometrical factor for meridional advection; geofac uses
!                   correct spherical geometry, and replaces acosp as the
!                   meridional geometrical factor in tpcore
!    dpi   : divergence at a grid point; used to calculate vertical motion (mb)
!    xmass : horizontal mass flux in E-W direction (mb)
!    ymass : horizontal mass flux in N-S direction (mb)
!
!  NOTES:
!
!******************************************************************************
!
      implicit none

      !----------------------
      !Argument declarations.
      !----------------------

      logical :: do_reduction
      real*8  :: geofac_pc
      real*8  :: geofac(ju1_gl:j2_gl)
      real*8  :: dpi   (i1:i2, ju1:j2, k1:k2)
      real*8  :: xmass (ilo:ihi, julo:jhi, k1:k2)
      real*8  :: ymass (ilo:ihi, julo:jhi, k1:k2)

      !----------------------
      !Variable declarations.
      !----------------------

      integer :: il, ij
      integer :: jst, jend

      !----------------
      !Begin execution.
      !----------------

      if (pr_diag) then
        Write (6,*) 'Calc_Divergence called by ', loc_proc
      end if

      !-------------------------
      !Calculate N-S divergence.
      !-------------------------

! No polar cap. (ccc, 8/3/10)
!      do ij = j1p, j2p
!
!        dpi(i1:i2,ij,:) =
!     &    (ymass(i1:i2,ij,:) - ymass(i1:i2,ij+1,:)) *
!     &    geofac(ij)
!
!      end do
!
!      if(j1p.ne.2) then
!        dpi(:,2,:) = 0.
!        dpi(:,j2-1,:) = 0.
!      endif

!      do ij = j1p_w, j2p_w
      do ij = j1p, j2p-1

        dpi(i1:i2,ij,:) =
     &    (ymass(i1:i2,ij,:) - ymass(i1:i2,ij+1,:)) *
     &    geofac(ij)

      end do

!-----lzh-----------------------
!      !===========================
!      call Do_Divergence_Pole_Sum
!      !===========================
!     &  (do_reduction, geofac_pc, dpi, ymass)
! comment out for nested NA (lzh, 07/20/2010)       
!        dpi(:,1,:) = 0.  ! (lzh, 07/20/2010)
!        dpi(:,j2,:) = 0.

      !-------------------------
      !Calculate E-W divergence.
      !-------------------------

      do ij = j1p,j2p
        do il = i1, i2-1
          dpi(il,ij,:) =
     &      dpi(il,ij,:) +
     &      xmass(il,ij,:) - xmass(il+1,ij,:)
        end do
! No periodicity. (ccc, 8/3/10)
!          dpi(i2,ij,:) =
!     &      dpi(i2,ij,:) +
!     &      xmass(i2,ij,:) - xmass(1,ij,:)
      end do

      ! Return to calling program
      END SUBROUTINE Calc_Divergence

!------------------------------------------------------------------------------

      SUBROUTINE Set_Press_Terms
     &  (dap, dbk, pres1, pres2, delp1, delpm, pu)
!
!******************************************************************************
! 
!  ROUTINE
!    Set_Press_Terms
!
!  AUTHORS
!    Philip Cameron-Smith and John Tannahill, GMI project @ LLNL (2003)
!
!  DESCRIPTION
!    This routine sets the pressure terms.
! 
!  ARGUMENTS
!    dap   : pressure difference across layer from (ai * pt) term (mb)
!    dbk   : difference in bi across layer - the dSigma term
!    pres1 : surface pressure at t1     (mb)
!    pres2 : surface pressure at t1+tdt (mb)
!    delp1 : pressure thickness, the psudo-density in a hydrostatic system
!            at t1 (mb)
!    delpm : pressure thickness, the psudo-density in a hydrostatic system
!            at t1+tdt/2 (approximate)  (mb)
!    pu    : pressure at edges in "u"   (mb)
!
!  NOTES:
!
!******************************************************************************
!
      implicit none
 
      !----------------------
      !Argument declarations.
      !----------------------

      real*8  :: dap  (k1:k2)
      real*8  :: dbk  (k1:k2)
      real*8  :: pres1(ilo:ihi, julo:jhi)
      real*8  :: pres2(ilo:ihi, julo:jhi)
      real*8  :: delp1(ilo:ihi, julo:jhi, k1:k2)
      real*8  :: delpm(ilo:ihi, julo:jhi, k1:k2)
      real*8  :: pu   (ilo:ihi, julo:jhi, k1:k2)

      !----------------------
      !Variable declarations.
      !----------------------

      integer :: il, ij, ik
      integer :: jst, jend

      !----------------
      !Begin execution.
      !----------------

      if (pr_diag) then
        Write (6,*) 'Set_Press_Terms called by ', loc_proc
      end if

      do ik = k1, k2

        delp1(:,:,ik) =
     &    dap(ik) + (dbk(ik) * pres1(:,:))

        delpm(:,:,ik) =
     &    dap(ik) + 
     &    (dbk(ik) * 0.5d0 * (pres1(:,:) + pres2(:,:)))

      end do

      do ij = ju1, j2
        do il = i1+1, i2
          pu(il,ij,:) =
     &      0.5d0 * (delpm(il,ij,:) + delpm(il-1,ij,:))
        end do

! No periodicity. (ccc, 8/3/10)
!          pu(i1,ij,:) =
!     &      0.5d0 * (delpm(i1,ij,:) + delpm(i2,ij,:))

      end do

      ! Return to calling program
      END SUBROUTINE Set_Press_Terms

!------------------------------------------------------------------------------

      SUBROUTINE XPAVG( P, IM )
!
!******************************************************************************
!  Subroutine XPAVG replaces each element of a vector with the average
!  of the entire array. (bmy, 5/7/03)
! 
!  Arguments as Input:
!  ============================================================================
!  (1 ) P  (REAL*8)  :: 1-D vector to be averaged
!  (2 ) IM (INTEGER) :: Dimension of P
!
!  Arguments as Output:
!  ============================================================================
!  (1 ) P  (REAL*8)  :: Contains average value of P in each element
!
!  NOTES:
!******************************************************************************
!
      ! References to F90 modules
      USE ERROR_MOD, ONLY : ERROR_STOP

      ! Arguments
      INTEGER, INTENT(IN)    :: IM
      REAL*8,  INTENT(INOUT) :: P(IM)

      ! Local variables
      REAL                   :: AVG

      !=================================================================
      ! XPAVG begins here!
      !=================================================================

      ! Error check IM
      IF ( IM == 0 ) THEN
         CALL ERROR_STOP( 'Div by zero!', 'XPAVG ("pjc_pfix_mod.f")' )
      ENDIF

      ! Take avg of entire P array
      AVG  = SUM( P ) / DBLE( IM ) 

      ! Store average value in all elements of P
      P(:) = AVG
 
      ! Return to calling program 
      END SUBROUTINE XPAVG

!------------------------------------------------------------------------------
      
      SUBROUTINE INIT_PJC_PFIX
!
!******************************************************************************
!  Subroutine INIT_PJC_PFIX allocates and initializes module arrays and
!  variables.  GMI dimension variables will be used for compatibility with
!  the Phil Cameron-Smith P-fixer. (bdf, bmy, 5/8/03)
!
!  NOTES:
!  01 Mar 2012 - R. Yantosca - Now use GET_AREA_M2(I,J,L) from grid_mod.F90
!  01 Mar 2012 - R. Yantosca - Now use GET_YMID_R(I,J,L) from grid_mod.F90
!******************************************************************************
!
      ! References to F90 modules
      USE GRID_MOD,     ONLY : GET_AREA_M2, GET_YMID_R
      USE ERROR_MOD,    ONLY : ALLOC_ERR,   ERROR_STOP
      USE PRESSURE_MOD, ONLY : GET_AP,      GET_BP

      USE CMN_SIZE_MOD  ! Size parameters
      USE CMN_GCTM_MOD  ! Re, PI, etc...

      ! Local variables
      INTEGER :: AS, I, J, L

      !=================================================================
      ! INIT_PJC_PFIX begins here!
      !
      ! Initialize dimensions for GMI pressure-fixer code 
      !=================================================================
      IMP_NBORDER = 0
      I1_GL       = 1 
      I2_GL       = IIPAR 
      JU1_GL      = 1 
      JV1_GL      = 1 
      J2_GL       = JJPAR 
      K1_GL       = 1 
      K2_GL       = LLPAR 
      ILO_GL      = I1_GL  - IMP_NBORDER 
      IHI_GL      = I2_GL  + IMP_NBORDER 
      JULO_GL     = JU1_GL - IMP_NBORDER 
      JVLO_GL     = JV1_GL - IMP_NBORDER 
      JHI_GL      = J2_GL  + IMP_NBORDER 
      I1          = I1_GL 
      I2          = I2_GL 
      JU1         = JU1_GL 
      JV1         = JV1_GL 
      J2          = J2_GL 
      K1          = K1_GL 
      K2          = K2_GL 
      ILO         = ILO_GL 
      IHI         = IHI_GL 
      JULO        = JULO_GL 
      JVLO        = JVLO_GL 
      JHI         = JHI_GL 
! No polar cap. (ccc, 8/3/10)
!      J1P         = 3 
      J1P         = 1 
      J2P         = J2_GL - J1P + 1 
! Used only to check dimensions
      ILAT        = J2_GL - JU1_GL + 1 
      ILONG       = I2_GL -  I1_GL + 1 
      IVERT       = K2_GL -  K1_GL + 1 

! To add a buffer zone to calculate p-fixer for nested grid 
! simulations. The p-fixer is not calculated for the edge boxes. 
! (lzh, ccc, 8/3/10)
      BUFF_SIZE     = 2
      I1_W          = I1_GL + BUFF_SIZE
      I2_W          = I2_GL - BUFF_SIZE
      JU1_W         = JU1_GL + BUFF_SIZE
      J2_W          = J2_GL - BUFF_SIZE 
      J1P_W         = 1 + BUFF_SIZE
      J2P_W         = J2_GL - J1P_W + 1

      ! Error check longitude
      IF ( ILONG /= IIPAR ) THEN
         CALL ERROR_STOP( 'Invalid longitude dimension ILONG!', 
     &                    'INIT_PJC_FIX ("pjc_pfix_mod.f")' )
      ENDIF
      
      ! Error check latitude
      IF ( ILAT /= JJPAR ) THEN
         CALL ERROR_STOP( 'Invalid latitude dimension ILAT!', 
     &                    'INIT_PJC_FIX ("pjc_pfix_mod.f")' )
      ENDIF

      ! Error check altitude
      IF ( IVERT /= LLPAR ) THEN
         CALL ERROR_STOP( 'Invalid altitude dimension IVERT!', 
     &                    'INIT_PJC_FIX ("pjc_pfix_mod.f")' )
      ENDIF
 
    
      !=================================================================      
      ! Allocate module arrays (use dimensions from GMI code)
      !=================================================================
      ALLOCATE( AI( K1_GL-1:K2_GL ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'AI' )

      ALLOCATE( BI( K1_GL-1:K2_GL ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'BI' )

      ALLOCATE( DAP( K1_GL:K2_GL ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'DAP' )

      ALLOCATE( DBK( K1_GL:K2_GL ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'DBK' )

      ALLOCATE( CLAT_FV( JU1_GL:J2_GL ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'CLAT_FV' )

      ALLOCATE( COSE_FV( JU1_GL:J2_GL+1 ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'COSE_FV' )

      ALLOCATE( COSP_FV( JU1_GL:J2_GL ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'COSP_FV' )

      ALLOCATE( DLAT_FV( JU1_GL:J2_GL ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'DLAT_FV' )

      ALLOCATE( ELAT_FV( JU1_GL:J2_GL+1 ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'ELAT_FV' )

      ALLOCATE( GEOFAC( JU1_GL:J2_GL ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'GEOFAC' )

      ALLOCATE( GW_FV( JU1_GL:J2_GL ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'GW_FV' )
      
      ALLOCATE( MCOR( I1_GL:I2_GL, JU1_GL:J2_GL ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'MCOR' )

      ALLOCATE( REL_AREA( I1_GL:I2_GL, JU1_GL:J2_GL ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'REL_AREA' )

      ALLOCATE( RGW_FV( JU1_GL:J2_GL ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'RGW_FV' )

      ALLOCATE( SINE_FV( JU1_GL:J2_GL+1 ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'SINE_FV' )


      !=================================================================
      ! Initialize arrays and variables
      !=================================================================

      ! Grid box surface areas [m2]
      DO J = JU1_GL, J2_GL
      DO I =  I1_GL, I2_GL
         MCOR(I,J) = GET_AREA_M2( I, J, 1 )
      ENDDO
      ENDDO

      ! Hybrid grid vertical coords: Ai [hPa] and Bi [unitless]
      DO L = K1_GL-1, K2_GL
         AI(L) = GET_AP( L+1 )
         BI(L) = GET_BP( L+1 )
      ENDDO

      ! Delta A [hPa] and Delta B [unitless]
      DO L = K1_GL, K2_GL
         !-------------------------------------------------------------
         ! NOTE:, this was the original code.  But since AI is already
         ! in hPa, we shouldn't need to multiply by PTOP again.  This
         ! should only matter for the fvDAS fields.  Also, DBK needs 
         ! to be positive (bmy, 5/8/03)
         !DAP(L) = ( AI(L) - AI(L-1) ) * PTOP   
         !DBK(L) = BI(L) - BI(L-1)
         !-------------------------------------------------------------
         DAP(L) = AI(L-1) - AI(L)  
         DBK(L) = BI(L-1) - BI(L)
      ENDDO
         
      ! Grid box center latitudes [radians]
      DO J = JU1_GL, J2_GL
         CLAT_FV(J) = GET_YMID_R( 1, J, 1 )
      ENDDO

      ! Longitude spacing
!clb  DLON_FV    = 2.d0 * PI / DBLE( 540 )        !(dan)
!clb
#if defined (GRID025x03125)
      DLON_FV    = 2.d0 * PI / DBLE( 1152 )
#else
      DLON_FV    = 2.d0 * PI / DBLE( 540 )
#endif
!clb
       
      ! Latitude edge at south pole [radians]
!      ELAT_FV(1) = -0.5d0 * PI 
! for nested NA or EA (lzh, 07/20/2010) 
      ELAT_FV(1) = CLAT_FV(1) - 0.25d0 * PI/DBLE(180)

      ! SIN and COS of lat edge at south pole [unitless]
!      SINE_FV(1) = -1.d0
!      COSE_FV(1) =  0.d0
! for nested NA or EA (lzh, 07/20/2010) 
      SINE_FV(1) =  SIN( ELAT_FV(1) )
      COSE_FV(1) =  COS( ELAT_FV(1) )
         
      ! Latitude edges [radians] (w/ SIN & COS) at intermediate latitudes
      DO J = JU1_GL+1, J2_GL  !2, JJPAR
         ELAT_FV(J) = 0.5d0 * ( CLAT_FV(J-1) + CLAT_FV(J) )
         SINE_FV(J) = SIN( ELAT_FV(J) )
         COSE_FV(J) = COS( ELAT_FV(J) )
      ENDDO

      ! Latitude edge at North Pole [radians]
!      ELAT_FV(J2_GL+1) = 0.5d0 * PI
! for nested NA or EA (lzh, 07/20/2010) 
      ELAT_FV(J2_GL+1) = CLAT_FV(J2_GL)+0.25d0* PI/DBLE(180)

      ! SIN of lat edge at North Pole
!      SINE_FV(J2_GL+1) = 1.d0
! for nested NA or EA (lzh, 07/20/2010) 
      SINE_FV(J2_GL+1) =  SIN( ELAT_FV(J2_GL+1) )
      COSE_FV(J2_GL+1) =  COS( ELAT_FV(J2_GL+1) )
       
      ! Latitude extent of South polar box [radians]
!      DLAT_FV(1) = 2.d0 * ( ELAT_FV(2) - ELAT_FV(1) ) 
! comment out for nested NA or EA (lzh, 07/20/2010)

      ! Latitude extent of boxes at intermediate latitudes [radians]
!      DO J = JU1_GL+1, J2_GL-1  ! 2, JJPAR-1
! for nested NA or EA (lzh, 07/20/2010) 
      DO J = JU1_GL, J2_GL
         DLAT_FV(J) = ELAT_FV(J+1) - ELAT_FV(J)
      ENDDO

      ! Latitude extent of North polar box [radians]
!      DLAT_FV(J2_GL) = 2.d0 * ( ELAT_FV(J2_GL+1) - ELAT_FV(J2_GL) ) 
! comment out for nested NA or EA (lzh, 07/20/2010)

      ! Other stuff
      DO J = JU1_GL, J2_GL
         GW_FV(J)   = SINE_FV(J+1) - SINE_FV(J)
         COSP_FV(J) = GW_FV(J)     / DLAT_FV(J)
         RGW_FV(J)  = 1.d0         / GW_FV(J)
      ENDDO

      ! Return to calling program
      END SUBROUTINE INIT_PJC_PFIX

!------------------------------------------------------------------------------

      SUBROUTINE CLEANUP_PJC_PFIX_GEOSFP_WINDOW
!
!******************************************************************************
!  Subroutine CLEANUP_PJC_PFIX deallocates all module arrays (bmy, 5/8/03)
!  
!  NOTES:
!  26 Sep 2013 - R. Yantosca - Renamed to CLEANUP_PJC_PFIX_GEOSFP_WINDOW
!******************************************************************************
!
      !=================================================================
      ! CLEANUP_PJC_PFIX begins here!
      !=================================================================
      IF ( ALLOCATED( AI       ) ) DEALLOCATE( AI       )
      IF ( ALLOCATED( BI       ) ) DEALLOCATE( BI       )
      IF ( ALLOCATED( CLAT_FV  ) ) DEALLOCATE( CLAT_FV  )
      IF ( ALLOCATED( COSE_FV  ) ) DEALLOCATE( COSE_FV  )
      IF ( ALLOCATED( COSP_FV  ) ) DEALLOCATE( COSP_FV  )
      IF ( ALLOCATED( DAP      ) ) DEALLOCATE( DAP      )
      IF ( ALLOCATED( DBK      ) ) DEALLOCATE( DBK      )
      IF ( ALLOCATED( DLAT_FV  ) ) DEALLOCATE( DLAT_FV  )
      IF ( ALLOCATED( ELAT_FV  ) ) DEALLOCATE( ELAT_FV  )
      IF ( ALLOCATED( GEOFAC   ) ) DEALLOCATE( GEOFAC   )
      IF ( ALLOCATED( GW_FV    ) ) DEALLOCATE( GW_FV    )
      IF ( ALLOCATED( MCOR     ) ) DEALLOCATE( MCOR     )
      IF ( ALLOCATED( REL_AREA ) ) DEALLOCATE( REL_AREA )
      IF ( ALLOCATED( RGW_FV   ) ) DEALLOCATE( RGW_FV   )
      IF ( ALLOCATED( SINE_FV  ) ) DEALLOCATE( SINE_FV )


      ! Return to calling program
      END SUBROUTINE CLEANUP_PJC_PFIX_GEOSFP_WINDOW

!------------------------------------------------------------------------------

      END MODULE PJC_PFIX_GEOSFP_WINDOW_MOD
